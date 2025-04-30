
library("rstudioapi", lib.loc="../../local/R/library/")
library("desc", lib.loc="../../local/R/library/")
library("withr", lib.loc="../../local/R/library/")
library("ps", lib.loc="../../local/R/library/")
library("usethis", lib.loc="../../local/R/library/")
library("devtools", lib.loc="../../local/R/library/")
library("tibble", lib.loc="../../local/R/library/")
library("readxl", lib.loc="../../local/R/library/")
library("scales", lib.loc="../../local/R/library/")

#' Only Keep one of AB, BA pairs from interactions
#'
#' @param data_input interaction matrix with gene names on row and column side
#' Assign the correct value for the pair in upper matrix, and assign the other one NA
AB_BA_removal <- function(data_input){
  
  # The interactions are not symmetric! (Only fill in the upper triangular matrix, and set the lower triangular matrix to NA)
  # (MC) For the yeast network, we filter out AB-BA pairs where the two interactions are both significant but have opposite sign. I suspect these are rare cases but should not be considered as real interactions. For those that are both significant and same sign, keep the pair with the more significant FDR.
  genes_AB_BA = sort(intersect(row.names(data_input), colnames(data_input))) # Common library and query genes
  
  if(length(genes_AB_BA) == 0)
    return(data_input)
  
  # For reference data (std/str cutoff), calcualte and keep value on upper triangular matrix only
  for (i in 1 : length(genes_AB_BA)){
    if (i == length(genes_AB_BA)){ next }
    
    for (j in (i+1) : length(genes_AB_BA)){
      
      # If both of them are zero, these will be skipped anyway
      if( (data_input[genes_AB_BA[i], genes_AB_BA[j]] == 0) & (data_input[genes_AB_BA[j], genes_AB_BA[i]] == 0) ){
        data_input[genes_AB_BA[j], genes_AB_BA[i]] <- NA
        next
      }
      
      # If one of them is zero, use the other
      if(data_input[genes_AB_BA[i], genes_AB_BA[j]] == 0){ 
        data_input[genes_AB_BA[i], genes_AB_BA[j]] <- data_input[genes_AB_BA[j], genes_AB_BA[i]]
        data_input[genes_AB_BA[j], genes_AB_BA[i]] <- NA
      } else{ # both have values
        
        # No need to change (already correct)
        if(data_input[genes_AB_BA[j], genes_AB_BA[i]] == 0){
          data_input[genes_AB_BA[j], genes_AB_BA[i]] <- NA
          next
        }
        
        # So, both have non-zero values
        # Different sign -> no interactions
        if(sign(data_input[genes_AB_BA[i], genes_AB_BA[j]]) != sign(data_input[genes_AB_BA[j], genes_AB_BA[i]])){
          data_input[genes_AB_BA[i], genes_AB_BA[j]] <- 0
          data_input[genes_AB_BA[j], genes_AB_BA[i]] <- NA
        } else{ # same sign -> Keep the larger (abs) one
          data_input[genes_AB_BA[i], genes_AB_BA[j]] <- sign(data_input[genes_AB_BA[i], genes_AB_BA[j]]) * max(abs(c(data_input[genes_AB_BA[i], genes_AB_BA[j]], data_input[genes_AB_BA[j], genes_AB_BA[i]])))
          data_input[genes_AB_BA[j], genes_AB_BA[i]] <- NA
        }
      }
    }
  }
  
  return(data_input)
}


###############################################
#========data load and preparation============#
###############################################

#source('Within_Between_Summary_or.R')

# load CORUM complex information
load("data_complex.rda")

# Read the qGI data (non-redundant)
#qGi scores
data_qGI <- read.table('qGI_20211111_nRed.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE) 
# FDR scores
data_fdr <- read.table('gi_FDR_20211111_nRed.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE) 

# Sort data rows by library gene names
names.genes <- row.names(data_qGI)
if (is.unsorted(names.genes)){
  ind <- order(names.genes)
  data_qGI <- data_qGI[ind,]
  data_fdr <- data_fdr[ind,]
}

# Sort data rows by query gene names
names.queries <- names(data_qGI)
if (is.unsorted(names.queries)){
  ind <- order(names.queries)
  data_qGI <- data_qGI[,ind]
  data_fdr <- data_fdr[,ind]
}

# Apply Standard cutoff (|qGI| >= 0.3 and FDR <= 0.1) to filter data
# Set interaction score to 0 if absolute qGI score is less than .3 or FDR is greater than .1
data_qGI_std <- data_qGI
data_qGI_std[(abs(data_qGI) < 0.3) | (data_fdr > 0.1)] <- 0 

# Stringent cutoff
# data_qGI_str <- data_qGI_nRed
# data_qGI_str[(abs(data_qGI_nRed) < 0.7) | (data_fdr_nRed > 0.05)] <- 0

# load core essential genes from Predicted_essential_genes.xlsx file, sheet 1 (Core_essential_genes)
data_input <- read_xlsx('Predicted_essential_genes.xlsx', sheet = 1, col_names=FALSE)
# get which essential genes are in the library gene set
genes_ess <- intersect(row.names(data_qGI_std), data_input$...1)

data_standard = data_complex # set data standard to CORUM complex
data_interaction = data_qGI_std # set interaction data to filtered qGI score data

# set minimum # of library and query gene required to be present in a complex to be selected to be in the analysis
min_no_of_genes_library = 3
min_no_of_genes_query = 1

postfix_out = 'Complex'

#remove TAPT1_375_rich to make the data non-redundant; there is another TAPT1 screen 
data_interaction <- data_interaction[,!names(data_interaction) %in% c("TAPT1_375_rich")]

				   
# Reduce the background to the set of domain (i.e. CORUM complex standard) genes
# get overlap of library genes and standard genes 
genes_in_standard <- intersect(row.names(data_interaction), unlist(strsplit(data_standard$Genes, split = ';')))
# filter interaction data to have only library genes that overlap with genes from the standard
data_interaction <- data_interaction[genes_in_standard, ]

#get the list of libary genes
genes_library <- row.names(data_interaction)

#get the list of query screens and extract the query gene names 
genes_queries_with_id <- names(data_interaction)
genes_queries <- unique(unlist(lapply(strsplit(genes_queries_with_id, '_'), '[[', 1) ))
## Some library and query names are different. Let's change the query names.
# KIAA0195 in the lib is called TMEM94 as a query and KIAA1432 in the lib is called RIC1 as a query.
# 'TMEM94' %in% genes_queries # TRUE
# 'RIC1' %in% genes_queries # TRUE
genes_queries[match(c('TMEM94', 'RIC1'), genes_queries)] <- c('KIAA0195', 'KIAA1432')
# reset column names to query gene names instead of screen names
names(data_interaction) <- genes_queries

############################################
#========background calculation============#
############################################

# Background / Population (full interaction space)
bg_num_actual_int_neg <- sum(data_interaction < 0, na.rm = T) # total number of negative GIs in the entire data
bg_num_actual_int_pos <- sum(data_interaction > 0, na.rm = T) # total number of positive GIs in the entire data
bg_num_actual_int_all <- bg_num_actual_int_pos + bg_num_actual_int_neg # total number of GIs (positive and negative) in the entire data
bg_num_tested_int <- dim(data_interaction)[1] * dim(data_interaction)[2] # total number of tested interactions

density_bg <- bg_num_actual_int_all / bg_num_tested_int # fraction of GIs in total tested interaction
density_neg_bg <- bg_num_actual_int_neg / bg_num_tested_int # fraction of negative GIs in total tested interaction
density_pos_bg <- bg_num_actual_int_pos / bg_num_tested_int # fraction of positive GIs in total tested interaction

#Convert self interactions (same library and query gene) to NA's
common_genes <- sort(intersect(genes_queries, row.names(data_interaction)))
for (gene in common_genes){
	data_interaction[gene, gene] <- NA
}

## ==================================================================================
##    Filter 1: Select complexes with a minimum number of library and query genes
## ==================================================================================

# for the modules (complex) in the standard, generate lists to store library and query gene overlaps and overlap size
# this is required to select complexes that meet the minimum requirement of library and query gene numbers

in_library <- numeric(length(data_standard$ID))
in_query <- numeric(length(data_standard$ID))
in_combined <- numeric(length(data_standard$ID))

library_list <- rep('', length(data_standard$ID))
query_list <- rep('', length(data_standard$ID))
unique_list <- rep('', length(data_standard$ID))

for (i in 1 : length(data_standard$ID)){
	gene_list = unlist(strsplit(data_standard$Genes[i], ';')) # To upper misses some of the genes?
	gene_list = gsub(' ', '', gene_list) # Replacing any spaces with empty string

	tmp_lib_genes <- intersect(gene_list, genes_library)  #overlap with library genes
	in_library[i] <- length(tmp_lib_genes) # size of overlap
	library_list[i] <- paste(tmp_lib_genes, collapse = ';') # Replacing any spaces with empty string

	tmp_query_genes <- intersect(gene_list, genes_queries) #overlap with query genes
	in_query[i] <- length(tmp_query_genes) # size of overlap
	query_list[i] <- paste(tmp_query_genes, collapse = ';') # Replacing any spaces with empty string

	tmp_combined <- union(tmp_lib_genes,tmp_query_genes) #overlap with library and query genes
	in_combined[i] <- length(tmp_combined) # size of overlap
	unique_list[i] <- paste(tmp_combined, collapse = ';') # Replacing any spaces with empty string
}

# generate a dataframe with the created lists
data_standard_lib_query <- cbind(data_standard, in_library, in_query, in_combined, query_genes = query_list, library_genes = library_list, library_query_genes = unique_list)

write.table(data_standard_lib_query, paste0('./output_library_3_or_query_1_complex/All_initial_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)

# ** Filter by number of genes on the library and query side
# select complex if it contains --
#	at least required minimum # of library genes(3) 
#	or at least required minimum # of query genes(1)  
#	or union of library and query gene overlap is at least required minimum # of library genes + required minimum # of query genes (4)
data_standard_lib_query <- data_standard_lib_query[((data_standard_lib_query$in_library >= min_no_of_genes_library) | (data_standard_lib_query$in_query >= min_no_of_genes_query) | data_standard_lib_query$in_combined >= (min_no_of_genes_library + min_no_of_genes_query)), ]

# Sort by number of genes on query side
ind <- order(data_standard_lib_query$in_query, decreasing = T)
data_standard_lib_query <- data_standard_lib_query[ind,] 

#generate dataframe with final list of complexes
data_module_cand <- data_standard_lib_query[,c('Name', 'Genes', 'Length', 'in_library', 'in_query', 'library_genes', 'query_genes')]
write.table(data_module_cand, paste0('./output_library_3_or_query_1_complex/All_Represented_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)

# remove duplicate rows (to remove same complexes that are annotated multiple times with different IDs)
data_module_cand <- unique(data_module_cand)
write.table(data_module_cand, paste0('./output_library_3_or_query_1_complex/All_Represented_unique_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)

## ======================================================================================
##        Generate data for within-between analysis
## ======================================================================================
# Size of dots = Density of interaction (# interaction / # possible pairs)
# (first gene from the library side, and second gene is from the query side)

#matrix to store # of tested interactions between modules (complexes)
module_no_tested_int <- matrix(NA, dim(data_module_cand)[1], dim(data_module_cand)[1])
  
#matrix to store interaction density between modules (complexes)
density_all <- matrix(NA, dim(data_module_cand)[1], dim(data_module_cand)[1])
row.names(density_all) <- data_module_cand$Name
colnames(density_all) <- data_module_cand$Name

#matrix to store interaction enrichment between modules (complexes) 
enrich_all <- matrix(NA, dim(data_module_cand)[1], dim(data_module_cand)[1])
row.names(enrich_all) <- data_module_cand$Name
colnames(enrich_all) <- data_module_cand$Name

#matrix to store # of actual interactions between modules (complexes) 
module_no_actual_int_all <- matrix(NA, dim(data_module_cand)[1], dim(data_module_cand)[1])
  
#matrix to store  Negative interactions - density, enrichment, # of observed negative GIs in modules
density_neg <- matrix(NA, dim(data_module_cand)[1], dim(data_module_cand)[1])
enrich_neg <- matrix(NA, dim(data_module_cand)[1], dim(data_module_cand)[1])
row.names(enrich_neg) <- data_module_cand$Name
colnames(enrich_neg) <- data_module_cand$Name
row.names(density_neg) <- data_module_cand$Name
colnames(density_neg) <- data_module_cand$Name  
module_no_actual_int_neg <- matrix(NA, dim(data_module_cand)[1], dim(data_module_cand)[1])
  
#matrix to store Positive interactions - density, enrichment, # of observed positive GIs in modules
density_pos <- matrix(NA, dim(data_module_cand)[1], dim(data_module_cand)[1])
enrich_pos <- matrix(NA, dim(data_module_cand)[1], dim(data_module_cand)[1])
row.names(enrich_pos) <- data_module_cand$Name
colnames(enrich_pos) <- data_module_cand$Name
row.names(density_pos) <- data_module_cand$Name
colnames(density_pos) <- data_module_cand$Name 
module_no_actual_int_pos <- matrix(NA, dim(data_module_cand)[1], dim(data_module_cand)[1])
  
# matrix to store total and fraction of essential genes in library module overlap
frac_ess_library <- matrix(NA, dim(data_module_cand)[1], dim(data_module_cand)[1])
ess_genes_library <- matrix(NA, dim(data_module_cand)[1], dim(data_module_cand)[1])

# create dataframe to store library gene, query gene, type of GI, # of GIs, module 1, module 2
all_pairs = data.frame(library = character(), query = character(), 
				   type = character(), GI = numeric(), 
				   domain1 = character(), domain2 = character(),
				   stringsAsFactors = F)

# for each selected module(complex)					   
for (i in 1 : dim(data_module_cand)[1]){
	print(i)
	# get list of module genes from ;-separated string
	gene_list_i = unlist(strsplit(data_module_cand$Genes[i], ';')) 
	# replace any whotespace with empty string
	gene_list_i = gsub(' ', '', gene_list_i) 
	
	# get library and query gene lists that are present in the current module (complex)
	# Current library genes: find overlap of libary and module genes and get indices of library gene list
	library_ind_i <- match(intersect(gene_list_i, genes_library), genes_library) 
	# Current Query genes: find overlap of query and module genes and get indices of query gene list
	query_ind_i <- match(intersect(gene_list_i, genes_queries), genes_queries) 
	
	## ======================================================================================
	##        Generate data for within complex analysis
	## ======================================================================================
	# For within analyses, values are stored across the diagonal of each matrix
	
	# if both library and query gene overlaps are greater than 0
	if(length(library_ind_i)>0 & length(query_ind_i)>0)
	{	
		# get subset of interaction matrix pertaining to current set of library and query genes
		Within_matrix <- data_interaction[library_ind_i, query_ind_i, drop = F]
		
		#store # of tested GIs in module: # of current library genes * # of current query genes
		module_no_tested_int[i,i]  <- dim(Within_matrix)[1] * dim(Within_matrix)[2]
		
		# AB-BA removal: Only Keep one of AB, BA pairs from interactions
		Within_matrix <- AB_BA_removal(Within_matrix) 
		# module_no_tested_int[i,i]  <- sum(!is.na(Within_matrix)) # Not counting self / NA
		
		# 0 and NA's are not counted
		#store # of observed negative GIs in module: # of GIs with negative scores in module
		module_no_actual_int_neg[i,i] <- sum(Within_matrix < 0, na.rm = T)
		#store density of negative GIs in module: # of observed negative GIs / # of total tested GIs
		density_neg[i,i] <- round(module_no_actual_int_neg[i,i] / module_no_tested_int[i,i], 4)
		# enrichment analysis  hyper-geometric test: test if # of negative GIs in module is significant
		enrich_neg[i,i]  <- phyper(module_no_actual_int_neg[i,i] - 1, bg_num_actual_int_neg, 
								   bg_num_tested_int - bg_num_actual_int_neg, module_no_tested_int[i,i], 
								   lower.tail = FALSE)
								   
		#store # of observed positive GIs in module: # of GIs with positive scores in module
		module_no_actual_int_pos[i,i] <- sum(Within_matrix > 0, na.rm = T)
		#store density of positive GIs in module: # of observed positive GIs / # of total tested GIs
		density_pos[i,i] <- round(module_no_actual_int_pos[i,i] / module_no_tested_int[i,i], 4)
		# enrichment analysis  hyper-geometric test: test if # of positive GIs in module is significant
		enrich_pos[i,i]  <- phyper(module_no_actual_int_pos[i,i] - 1, bg_num_actual_int_pos, 
								   bg_num_tested_int - bg_num_actual_int_pos, module_no_tested_int[i,i], 
								   lower.tail = FALSE)
		
		#store # of observed GIs in module: # of GIs with positive scores in module + # of GIs with negative scores in module
		module_no_actual_int_all[i,i] <- module_no_actual_int_neg[i,i] + module_no_actual_int_pos[i,i]
		#store density of GIs in module: # of observed GIs / # of total tested GIs
		density_all[i,i] <- round(module_no_actual_int_all[i,i] / module_no_tested_int[i,i], 4)
		# enrichment analysis  hyper-geometric test: test if # of GIs in module is significant
		enrich_all[i,i]  <- phyper(module_no_actual_int_all[i,i] - 1, bg_num_actual_int_all, 
								   bg_num_tested_int - bg_num_actual_int_all, module_no_tested_int[i,i], 
								   lower.tail = FALSE)
		
		# Store fraction of essential genes in library gene set in module
		overlap_library <- genes_library[library_ind_i]		
		frac_ess_library[i,i] <- round(length(intersect(overlap_library, genes_ess)) / length(overlap_library), 2)
		ess_genes_library[i,i] <- paste(intersect(overlap_library, genes_ess), collapse = ';')
		
		# Store all neg/pos gene pairs contributing to within analysis
		# store library gene, query gene, type of GI, # of GIs, module 1, module 2 in dataframe
		log_ind <- !is.na(Within_matrix) & Within_matrix != 0
		for (row_ind in 1:dim(log_ind)[1]){
		  for (col_ind in 1:dim(log_ind)[2]){
			if(log_ind[row_ind,col_ind]){
				all_pairs = rbind(all_pairs, data.frame(
				library = row.names(Within_matrix)[row_ind], 
				query = colnames(Within_matrix)[col_ind], 
				type = 'Within', 
				GI = round(Within_matrix[row_ind, col_ind],2), 
				domain1 = data_module_cand$Name[i], 
				domain2 = '',
				stringsAsFactors = F))
			}
		  }
		}
	}

	# if processing the last module, skip the between analysis. This module has already been compared to all other modules for between module analysis 
	if(i == dim(data_module_cand)[1]) next
	
	## ======================================================================================
	##        Generate data for between module (complex) analysis
	## ======================================================================================
	# There maybe some between pairs that are also within pairs (for either complex i or j) and we have to remove those pairs from between pair analysis
	
	# for modules listed after the current module
	for (j in (i+1) : dim(data_module_cand)[1]){
		# get list of module j genes from ;-separated string
	  gene_list_j = unlist(strsplit(data_module_cand$Genes[j], ';'))
	  # replace any whotespace with empty string
	  gene_list_j = gsub(' ', '', gene_list_j)
	  
	  # Unique library genes and unique query genes from module j
	  library_ind_j <- match(intersect(gene_list_j, genes_library), genes_library)
	  query_ind_j <- match(intersect(gene_list_j, genes_queries), genes_queries)
	  
	  # Within + between space
	  # get subset of interaction matrix: rows (library genes from module i and j), columns (query genes from module i and j)
	  WB_matrix <- data_interaction[union(library_ind_i, library_ind_j), union(query_ind_i, query_ind_j), drop = F]
	  WB_matrix[genes_library[library_ind_i], genes_queries[query_ind_i]] <- NA # Within GIs for module i set to NA
	  WB_matrix[genes_library[library_ind_j], genes_queries[query_ind_j]] <- NA # Within GIs for module j set to NA
	  # WB_matrix <- AB_BA_removal(WB_matrix) # Don't need this as AB/BA are relevant for within only
	  
	  # Removing within-pairs (by removing common genes from library side)
	  # library_ind_i_not_j <- setdiff(library_ind_i, library_ind_j)
	  # library_ind_j_not_i <- setdiff(library_ind_j, library_ind_i)
	  # query_ind_i_not_j <- setdiff(query_ind_i, query_ind_j)
	  # query_ind_j_not_i <- setdiff(query_ind_j, query_ind_i)
	  # module_no_tested_int[i,j] <- length(library_ind_i_not_j) * length(query_ind_j_not_i) + length(library_ind_j_not_i) * length(query_ind_i_not_j)
	  
	  #store # of tested GIs in module (i,j) pair: all GIs except self-module interactions (Self-interactions should already be considered in within-module analysis)
	  module_no_tested_int[i,j] <- sum(!is.na(WB_matrix)) 
	  module_no_tested_int[j,i] <- module_no_tested_int[i,j]
	  
	  # If the candidate interactions list is not empty
	  if (module_no_tested_int[i,j] > 0){
		
		## ---------- Way3: Symmetric ----------
		# (Library side i * Query side j + Library side j * Query side i)
		# Self interaction are usually NA here
		# Common interactions between Modules are not removed
		
		# === neg ===
		# module_no_actual_int_neg[i,j] <- sum(data_interaction[library_ind_j_not_i, query_ind_i_not_j] < 0, na.rm = T) + sum(data_interaction[library_ind_i_not_j, query_ind_j_not_i] < 0, na.rm = T)
		
		#store # of observed negative GIs in module (i,j) pair: # of GIs with negative scores
		module_no_actual_int_neg[i,j] <- sum(WB_matrix < 0, na.rm = TRUE)
		#store density of negative GIs in module (i,j) pair: # of observed negative GIs / # of total tested GIs
		density_neg[i,j] <- round(module_no_actual_int_neg[i,j] / module_no_tested_int[i,j], 4)
		# enrichment analysis  hyper-geometric test: test if # of negative GIs in  module (i,j) pair is significant
		enrich_neg[i,j]  <- phyper(module_no_actual_int_neg[i,j] - 1, bg_num_actual_int_neg,
								   bg_num_tested_int - bg_num_actual_int_neg, module_no_tested_int[i,j],
								   lower.tail = FALSE)
		module_no_actual_int_neg[j,i] <- module_no_actual_int_neg[i,j]
		density_neg[j,i] <- density_neg[i,j]
		enrich_neg[j,i] <- enrich_neg[i,j]
		
		# === pos ===
		# module_no_actual_int_pos[i,j] <- sum(data_interaction[library_ind_j_not_i, query_ind_i_not_j] > 0, na.rm = T) + sum(data_interaction[library_ind_i_not_j, query_ind_j_not_i] > 0, na.rm = T)
		
		#store # of observed positive GIs in module (i,j) pair: # of GIs with positive scores
		module_no_actual_int_pos[i,j] <- sum(WB_matrix > 0, na.rm = TRUE)
		#store density of positive GIs in module (i,j) pair: # of observed positive GIs / # of total tested GIs
		density_pos[i,j] <- round(module_no_actual_int_pos[i,j] / module_no_tested_int[i,j], 4)
		# enrichment analysis  hyper-geometric test: test if # of positive GIs in  module (i,j) pair is significant
		enrich_pos[i,j]  <- phyper(module_no_actual_int_pos[i,j] - 1, bg_num_actual_int_pos, 
								   bg_num_tested_int - bg_num_actual_int_pos, module_no_tested_int[i,j], 
								   lower.tail = FALSE) 		
		module_no_actual_int_pos[j,i] <- module_no_actual_int_pos[i,j]
		density_pos[j,i] <- density_pos[i,j]
		enrich_pos[j,i] <- enrich_pos[i,j]
		
		#store # of observed GIs in module (i,j) pair: # of GIs with positive scores + # of GIs with negative scores
		module_no_actual_int_all[i,j] <- module_no_actual_int_neg[i,j] + module_no_actual_int_pos[i,j]
		#store density of GIs in module (i,j) pair: # of observed GIs / # of total tested GIs
		density_all[i,j] <- round(module_no_actual_int_all[i,j] / module_no_tested_int[i,j], 4)
		# enrichment analysis  hyper-geometric test: test if # of GIs in  module (i,j) pair is significant
		enrich_all[i,j]  <- phyper(module_no_actual_int_all[i,j] - 1, bg_num_actual_int_all, 
								   bg_num_tested_int - bg_num_actual_int_all, module_no_tested_int[i,j], 
								   lower.tail = FALSE)		
		module_no_actual_int_all[j,i] <- module_no_actual_int_all[i,j]
		density_all[j,i] <- density_all[i,j]
		enrich_all[j,i] <- enrich_all[i,j]
		
		# Store fraction of essential genes in library gene set in module
		overlap_library <- genes_library[union(library_ind_i, library_ind_j)]		
		frac_ess_library[i,j] <- round(length(intersect(overlap_library, genes_ess)) / length(overlap_library), 2)
		ess_genes_library[i,j] <- paste(intersect(overlap_library, genes_ess), collapse = ';')
		frac_ess_library[j,i] <- frac_ess_library[i,j]
		ess_genes_library[j,i] <- ess_genes_library[i,j]
		
		#  Store all neg/pos gene pairs contributing to between analysis
		# store library gene, query gene, type of GI, # of GIs, module 1, module 2 in dataframe
		log_ind <- !is.na(WB_matrix) & WB_matrix != 0
		for (row_ind in 1:dim(log_ind)[1]){
		  for (col_ind in 1:dim(log_ind)[2]){
			if(log_ind[row_ind,col_ind]){
			  all_pairs = rbind(all_pairs, data.frame(
				library = row.names(WB_matrix)[row_ind], 
				query = colnames(WB_matrix)[col_ind], 
				type = 'Between', 
				GI = round(WB_matrix[row_ind, col_ind],2), 
				domain1 = data_module_cand$Name[i], 
				domain2 = data_module_cand$Name[j],
				stringsAsFactors = F))
			}
		  }
		}		
	  }
	}
}
  
# write to file: matrix with gene pairs contributing to analysis
write.table(all_pairs, './output_library_3_or_query_1_complex/All_WB_interaction_pairs_3_or_1.txt', sep = '\t', col.names = TRUE, row.names = FALSE)
  
# write to file: matrix with module enrichment scores
write.table(enrich_neg, paste0('./output_library_3_or_query_1_complex/enrich_neg_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)
write.table(enrich_pos, paste0('./output_library_3_or_query_1_complex/enrich_pos_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)
write.table(enrich_all, paste0('./output_library_3_or_query_1_complex/enrich_all_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)

# write to file: matrix with module interaction counts
write.table(data_module_cand, paste0('./output_library_3_or_query_1_complex/data_module_cand_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)
write.table(module_no_tested_int, paste0('./output_library_3_or_query_1_complex/module_no_tested_int_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)
write.table(module_no_actual_int_all, paste0('./output_library_3_or_query_1_complex/module_no_actual_int_all_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)
write.table(module_no_actual_int_neg, paste0('./output_library_3_or_query_1_complex/module_no_actual_int_neg_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)
write.table(module_no_actual_int_pos, paste0('./output_library_3_or_query_1_complex/module_no_actual_int_pos_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)

# write to file: matrix with module density		   
write.table(density_neg, paste0('./output_library_3_or_query_1_complex/density_neg_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)
write.table(density_pos, paste0('./output_library_3_or_query_1_complex/density_pos_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)
write.table(density_all, paste0('./output_library_3_or_query_1_complex/density_all_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)

# write to file: matrix with essential gene fraction
write.table(frac_ess_library, paste0('./output_library_3_or_query_1_complex/frac_ess_library_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)
write.table(ess_genes_library, paste0('./output_library_3_or_query_1_complex/ess_genes_library_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)

## ======================================================================================
##        Compile and post-process within-between analysis data
## ======================================================================================

## ================ Adjust pvalues using benjamini hochberg method (Separately for within and between) ================
# Adjust enrichment analyis (hyper-geometric test) p-values from all GIs
diag(enrich_all) <- p.adjust(diag(enrich_all), "BH") # Within
enrich_all[upper.tri(enrich_all)] <- p.adjust(enrich_all[upper.tri(enrich_all)], "BH") # Between (upper matrix)
enrich_all[lower.tri(enrich_all)] <- p.adjust(enrich_all[lower.tri(enrich_all)], "BH") # Between (lower matrix)
# Adjust enrichment analyis (hyper-geometric test) p-values from negative GIs
diag(enrich_neg) <- p.adjust(diag(enrich_neg), "BH") # Within
enrich_neg[upper.tri(enrich_neg)] <- p.adjust(enrich_neg[upper.tri(enrich_neg)], "BH") # Between (upper matrix)
enrich_neg[lower.tri(enrich_neg)] <- p.adjust(enrich_neg[lower.tri(enrich_neg)], "BH") # Between (lower matrix)
# Adjust enrichment analyis (hyper-geometric test) p-values from positive GIs
diag(enrich_pos) <- p.adjust(diag(enrich_pos), "BH") # Within
enrich_pos[upper.tri(enrich_pos)] <- p.adjust(enrich_pos[upper.tri(enrich_pos)], "BH") # Between (upper matrix)
enrich_pos[lower.tri(enrich_pos)] <- p.adjust(enrich_pos[lower.tri(enrich_pos)], "BH") # Between (upper matrix) 
  
## ================ Write the output in a long format ================
#  Format for the Supplementary file for this enrichment
#  Row_Complex Column_Complex 
#  Type (within/between) 
#  No_of_interaction_tested	No_total_interactions	No_Negative_interactions	No_Positive_interactions 
#  Fraction_neg_to_pos	Purity_score 
#  Enrichment_pvalue_total	Enrichment_pvalue_Negative	Enrichment_pvalue_Positive		
#  density_neg	density_pos	density_all
#  fraction_of_essential_in_library	essential_genes (only for within)

#  Data for dot plot as a list of matrices: density and enrichment

#  create dataframe place-holder
data_output <- data.frame(Row_Complex = character(), 
						Column_Complex = character(),
						type = character(), 
						No_of_interaction_tested = numeric(), 
						No_total_interactions = numeric(),
						No_neg_interactions = numeric(),
						No_pos_interactions = numeric(),
						Fraction_neg_to_pos = numeric(),
						Purity_score = numeric(),
						Enrichment_total = numeric(),
						Enrichment_neg = numeric(),
						Enrichment_pos = numeric(),
						density_neg = numeric(),
						density_pos = numeric(),
						density_all = numeric(),
						fraction_of_essential_in_library = numeric(),
						essential_genes = character())
						
# populate dataframe from previously generated matrices
for (i in 1 : dim(data_module_cand)[1]){
	print(i)
	# We need this to print the last within too
	for (j in i : dim(data_module_cand)[1]){	  
	  # We do not have candidate here anyway, so skip!
	  # Everything passes?
	  if(is.na(module_no_tested_int[i,j])) next
	  if(module_no_tested_int[i,j] == 0) next
	  
	  data_output <- rbind(data_output, data.frame(Row_Complex = data_module_cand$Name[i], # module 1
												   Column_Complex = data_module_cand$Name[j], # module 2
												   type = if (i == j) 'within' else 'between', # type of analysis
												   No_of_interaction_tested = module_no_tested_int[i,j], # total # of tested GIs
												   No_total_interactions = module_no_actual_int_all[i,j], # total # of observed GIs
												   No_neg_interactions = module_no_actual_int_neg[i,j], # total # of observed negative GIs
												   No_pos_interactions = module_no_actual_int_pos[i,j], # total # of observed positieve GIs
												   Fraction_neg_to_pos = round(module_no_actual_int_neg[i,j] / module_no_actual_int_pos[i,j], 2), # (total # of observed negative GIs / total # of observed positieve GIs) 
												   Purity_score = round(module_no_actual_int_pos[i,j] / module_no_actual_int_all[i,j], 2), # (total # of observed positive GIs / total # of observed GIs)
												   Enrichment_total = enrich_all[i,j], # enrichment of module(s) for GIs
												   Enrichment_neg = enrich_neg[i,j], # enrichment of module(s) for negative GIs 
												   Enrichment_pos = enrich_pos[i,j], # enrichment of module(s) for positive GIs 
												   density_neg = density_neg[i,j], # density of positive GIs in module
												   density_pos = density_pos[i,j], # density of negative GIs in module
												   density_all = density_all[i,j], # density of GIs in module
												   fraction_of_essential_in_library = frac_ess_library[i,j], # fraction of essential genes in module (within only)
												   ess_genes_library = ess_genes_library[i,j])) # essential genes in module (within only)
	}
}

# write to file
write.table(data_output, file = paste0('./output_library_3_or_query_1_complex/Within_Between_Enrichment_', postfix_out, '.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)

## ================ Rescale purity scores ================
# replace infinite/nan purity scores with NA
data_output[is.infinite(data_output$Purity_score), c('Purity_score')] <- NA
data_output[is.nan(data_output$Purity_score), c('Purity_score')] <- NA
# rescale purity score to range between [-1,1]: -1 all negative GIs, 1 all positive GIs
data_output$Purity_score <- rescale(data_output$Purity_score, to = c(-1, 1))

# Convert factors to characters
i <- sapply(data_output, is.factor)
data_output[i] <- lapply(data_output[i], as.character)

# Make the significant enrichment driven by 1 pair insignificant (FDR of 1)
data_output[data_output$No_total_interactions < 2, ]$Enrichment_total = 1
data_output[data_output$No_neg_interactions < 2, ]$Enrichment_neg = 1
data_output[data_output$No_pos_interactions < 2, ]$Enrichment_pos = 1

# Sort by enrichment (total)
ind_enrichment <- order(data_output$Enrichment_total, decreasing = F)
data_output <- data_output[ind_enrichment,]  

# Sort by number of total interactions (decreasing)
ind_total <- order(data_output$No_total_interactions, decreasing = T)
data_output <- data_output[ind_total,]

# Sort by differences in number of interactions (negative vs positive)
ind_diff <- order(abs(data_output$No_neg_interactions - data_output$No_pos_interactions), decreasing = T)
data_output <- data_output[ind_diff,]

# Sort by within and then between
ind_type <- order(data_output$type, decreasing = TRUE) # We want all 'within' first, followed by 'between'
data_output <- data_output[ind_type,]

#Write to file
write.table(data_output, file = paste0('./output_library_3_or_query_1_complex/Within_Between_Enrichment_adjusted_', postfix_out, '_no_filter_purity_modified.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)

  