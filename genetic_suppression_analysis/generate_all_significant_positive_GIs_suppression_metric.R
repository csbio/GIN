# author: Arshia Z. Hassan. hassa418@umn.edu

local_R_path = "../../../local/R/library/"
library("reshape2", lib.loc=local_R_path)
library("tibble", lib.loc=local_R_path)
library("readxl", lib.loc=local_R_path)
library("dplyr", lib.loc=local_R_path)

# load qGI score and FDR data - non-rudandant screens
data_qGI_nRed <- read.table('./qGI_20211111_nRed.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE) 
data_fdr_nRed <- read.table('./gi_FDR_20211111_nRed.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE) 

data_qGI_nRed_melt <- setNames(melt(as.matrix(data_qGI_nRed)), c('library_gene', ' query_name', 'qGI_score'))
data_fdr_nRed_melt <- setNames(melt(as.matrix(data_fdr_nRed)), c('library_gene', ' query_name', 'FDR'))

all_significant_positive_GIs = merge(data_qGI_nRed_melt,data_fdr_nRed_melt,by = c('library_gene', ' query_name'))

split_screen_name<-function(x)
{
  x1<-strsplit(x,"_")[[1]][1] #split x with delimiter \. and retrieve first sub-string
  return( x1)
}

#setNames(melt(data_qGI_nRed), c('library_gene', ' query_name', 'qGI_score'))

# get all significant positive GIs with FDR < .1 and qGI score >.3
all_significant_positive_GIs = all_significant_positive_GIs[which(all_significant_positive_GIs$FDR<.1 & all_significant_positive_GIs$qGI_score>.3),]
all_significant_positive_GIs$query_gene <- unlist(lapply(as.character(all_significant_positive_GIs$` query_name`),split_screen_name) ) # extract query gene name from screen name

# read gene information 
centroid_file = "./Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx"
supp_table_2 = read_excel(centroid_file)

# add columns to all significant positive GIs indicating if the library gene and the query genes are hap1 essentials
all_significant_positive_GIs$is_library_HAP1_essential = all_significant_positive_GIs$library_gene %in% supp_table_2$`Gene Symbol`[which(supp_table_2$`HAP1_essential gene`==1)]
all_significant_positive_GIs$is_query_HAP1_essential = all_significant_positive_GIs$query_gene %in% supp_table_2$`Gene Symbol`[which(supp_table_2$`HAP1_essential gene`==1)]

# for each significant positive GIs, get single mutant fitness scores of library and query genes for both minimal and rich medias
all_significant_positive_GIs$library_SMF_minimal = supp_table_2$`HAP1_Mutant fitness_minimal (log2 fold change)`[match(all_significant_positive_GIs$library_gene,supp_table_2$`Gene Symbol`)]
all_significant_positive_GIs$library_SMF_rich = supp_table_2$`HAP1_Mutant fitness_rich (log2 fold change)`[match(all_significant_positive_GIs$library_gene,supp_table_2$`Gene Symbol`)]
all_significant_positive_GIs$query_SMF_minimal = supp_table_2$`HAP1_Mutant fitness_minimal (log2 fold change)`[match(all_significant_positive_GIs$query_gene,supp_table_2$`Gene Symbol`)]
all_significant_positive_GIs$query_SMF_rich = supp_table_2$`HAP1_Mutant fitness_rich (log2 fold change)`[match(all_significant_positive_GIs$query_gene,supp_table_2$`Gene Symbol`)]

write.csv(all_significant_positive_GIs,"./all_significant_positive_GIs.csv",row.names=F)

# remove rows from screen 'TAPT1_375_rich' and self intearction GIs (same library and query gene)
all_significant_positive_GIs = all_significant_positive_GIs[-which(all_significant_positive_GIs$` query_name`=='TAPT1_375_rich'),]
all_significant_positive_GIs = all_significant_positive_GIs[-which(all_significant_positive_GIs$library_gene==all_significant_positive_GIs$query_gene),]

get_screen_type<-function(x)
{
  x1<-strsplit(x,"_")[[1]][3] #split x with delimiter \. and retrieve first sub-string
  return( x1)
}
all_significant_positive_GIs$query_type = unlist(lapply(as.character(all_significant_positive_GIs$` query_name`),get_screen_type) )

# get list of mitochondrial genes
mito_gene_list = read.csv('./Mitochondial_genelist_1_26_2021_genes.tsv',header=F)

# add columns to all significant positive GIs indicating if the library and the query genes are mitochondrial genes
all_significant_positive_GIs$is_library_mito = rep(c(FALSE),times=nrow(all_significant_positive_GIs))
all_significant_positive_GIs$is_query_mito = rep(c(FALSE),times=nrow(all_significant_positive_GIs))
for( i in 1:nrow(all_significant_positive_GIs))
{
  if(all_significant_positive_GIs[i,"library_gene"] %in% mito_gene_list$V1)
  {
    all_significant_positive_GIs[i,'is_library_mito'] = T
  }
  if(all_significant_positive_GIs[i,"query_gene"] %in% mito_gene_list$V1)
  {
    all_significant_positive_GIs[i,'is_query_mito'] = T
  }
  
}

# set NA values to 0 for single mutant fitness scores
all_significant_positive_GIs$library_SMF_minimal[is.na(all_significant_positive_GIs$library_SMF_minimal)] <- 0
all_significant_positive_GIs$query_SMF_minimal[is.na(all_significant_positive_GIs$query_SMF_minimal)] <- 0
all_significant_positive_GIs$library_SMF_rich[is.na(all_significant_positive_GIs$library_SMF_rich)] <- 0
all_significant_positive_GIs$query_SMF_rich[is.na(all_significant_positive_GIs$query_SMF_rich)] <- 0

# for each significant positive GI, calculate minimum of the library and query genes for both rich and minimal media type
all_significant_positive_GIs = all_significant_positive_GIs %>% 
  rowwise() %>%
  mutate(minimal_SMF_min =min(library_SMF_minimal,query_SMF_minimal))

all_significant_positive_GIs = all_significant_positive_GIs %>% 
  rowwise() %>%
  mutate(rich_SMF_min =min(library_SMF_rich,query_SMF_rich))

# for each significant positive GI, calculate (library SMF + query SMF + qGI score) for both rich and minimal media type

all_significant_positive_GIs$minimal_SMF_QGI = all_significant_positive_GIs$library_SMF_minimal + 
  all_significant_positive_GIs$query_SMF_minimal + 
  all_significant_positive_GIs$qGI_score

all_significant_positive_GIs$rich_SMF_QGI = all_significant_positive_GIs$library_SMF_rich + 
  all_significant_positive_GIs$query_SMF_rich + 
  all_significant_positive_GIs$qGI_score

# for each significant positive GI, calculate ((library SMF + query SMF + qGI score) - minimum SMF / absolute(minimum SMF)) 
# based on the rich and minimal media type of the query gene
# given minimum SMF < -0.9 and minimum SMF < (library SMF + query SMF + qGI score)
all_significant_positive_GIs$score = 0
SMF_min_cutoff = -0.9 #-1
for(i in 1:nrow(all_significant_positive_GIs)){
  if(all_significant_positive_GIs[i,"query_type"]=='min')
  {
    SMF_QGI = all_significant_positive_GIs[i,"minimal_SMF_QGI"]
    SMF_min = all_significant_positive_GIs[i,"minimal_SMF_min"]
  }
  else if(all_significant_positive_GIs[i,"query_type"]=='rich')
  {
    SMF_QGI = all_significant_positive_GIs[i,"rich_SMF_QGI"]
    SMF_min = all_significant_positive_GIs[i,"rich_SMF_min"]
  }
  if(SMF_min<(SMF_min_cutoff) & SMF_min<SMF_QGI){
    all_significant_positive_GIs[i,'score'] = ((SMF_QGI - SMF_min)/abs(SMF_min))
  }
}

write.csv(all_significant_positive_GIs,"./all_significant_positive_GIs_supp_metric_2_.9_.csv",row.names=F)

