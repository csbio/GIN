# author: Arshia Z. Hassan. hassa418@umn.edu

local_R_path = "../../local/R/library/"
library("tidyr", lib.loc=local_R_path)
library("stringr", lib.loc=local_R_path)
library("dplyr", lib.loc=local_R_path)
library("Matrix", lib.loc=local_R_path)
library("reshape2", lib.loc=local_R_path)
library("readxl", lib.loc=local_R_path)
library("writexl", lib.loc=local_R_path)
library("withr", lib.loc=local_R_path)
library("RColorBrewer", lib.loc=local_R_path)
library("labeling", lib.loc=local_R_path)
library("farver", lib.loc=local_R_path)
library("digest", lib.loc=local_R_path)
library("ggplot2", lib.loc=local_R_path)
library("ggthemes", lib.loc=local_R_path)
library("gridExtra", lib.loc=local_R_path)

# load OMOIM data 
OMIM_folder = "./"
morbidmap_file = "morbidmap_.txt"
morbidmap = read.csv(file.path(OMIM_folder, morbidmap_file), sep = '\t')
# extract disease gene list
morbidmap_ = morbidmap %>%
  mutate(Gene.Locus.And.Other.Related.Symbols = strsplit(as.character(Gene.Locus.And.Other.Related.Symbols), ", ")) %>%
  unnest(Gene.Locus.And.Other.Related.Symbols)

# get non-essential gene list with mutant fitness log2 fold change score greater than -.5
supp_table_2_file = "./Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx"
supp_table_2 = read_excel(supp_table_2_file)
get_screen_type<-function(x)
{
  x1<-strsplit(x,"_")[[1]][3] #split x with delimiter \. and retrieve first sub-string
  return( x1)
}
NE_genes = supp_table_2$`Gene Symbol`[which(supp_table_2$`HAP1_Mutant fitness_rich (log2 fold change)`> (-.5) & is.na(supp_table_2$`HAP1_essential gene`)==T)]

# get non-essential OMIM genes 
NE_OMIM_genes = unique(intersect(NE_genes,morbidmap_$Gene.Locus.And.Other.Related.Symbols))
write.csv(NE_OMIM_genes,"./output/NE_OMIM_genes.csv",row.names=F)

# generate Synthetic Lethal data

# load qGI score and FDR data - non-rudandant screens
data_qGI_nRed <- read.table('./qGI_20211111_nRed.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
data_fdr_nRed <- read.table('./gi_FDR_20211111_nRed.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# get significant negative GIs with FDR<.1 and qGI score < 0
data_qGI_nRed_melt <- setNames(melt(as.matrix(data_qGI_nRed)), c('library_gene', ' query_name', 'qGI_score'))
data_fdr_nRed_melt <- setNames(melt(as.matrix(data_fdr_nRed)), c('library_gene', ' query_name', 'FDR'))

all_significant_neg_GIs = merge(data_qGI_nRed_melt,data_fdr_nRed_melt,by = c('library_gene', ' query_name'))

# remove rows with query screen 'TAPT1_375_rich'
all_significant_neg_GIs = all_significant_neg_GIs[-which(all_significant_neg_GIs$` query_name`=='TAPT1_375_rich'),] 

all_significant_neg_GIs = all_significant_neg_GIs[which(all_significant_neg_GIs$FDR<0.1 & all_significant_neg_GIs$qGI_score<0),]

# keep significant negative GIs where library gene in not hap1-essential
all_significant_neg_GIs = all_significant_neg_GIs[which(all_significant_neg_GIs$library_gene %in% NE_OMIM_genes),]

# extract query gene name from query screen name
split_screen_name<-function(x)
{
  x1<-strsplit(x,"_")[[1]][1] #split x with delimiter \. and retrieve first sub-string
  return( x1)
}
all_significant_neg_GIs$query_gene <- unlist(lapply(as.character(all_significant_neg_GIs$` query_name`),split_screen_name) )

# remove self-interactions from significant negative GIs
all_significant_neg_GIs = all_significant_neg_GIs[-which(all_significant_neg_GIs$query_gene==all_significant_neg_GIs$library_gene),]

# extract query type from query screen name
get_screen_type<-function(x)
{
  x1<-strsplit(x,"_")[[1]][3] #split x with delimiter \. and retrieve first sub-string
  return( x1)
}
all_significant_neg_GIs$query_type = unlist(lapply(as.character(all_significant_neg_GIs$` query_name`),get_screen_type) )

# get single mutant fitness score based on query type (rich or minimal)
all_significant_neg_GIs$SMF_library = 0
all_significant_neg_GIs$SMF_query = 0
for(i in 1:nrow(all_significant_neg_GIs))
{

  if(all_significant_neg_GIs[i,'query_type']=='rich')
  {
    all_significant_neg_GIs[i,'SMF_library'] = supp_table_2$`HAP1_Mutant fitness_rich (log2 fold change)`[which(supp_table_2$`Gene Symbol`== all_significant_neg_GIs[i,"library_gene"])]
    if(all_significant_neg_GIs[i,"query_gene"] %in% supp_table_2$`Gene Symbol`)
    {
      all_significant_neg_GIs[i,'SMF_query'] = supp_table_2$`HAP1_Mutant fitness_rich (log2 fold change)`[which(supp_table_2$`Gene Symbol`== all_significant_neg_GIs[i,"query_gene"])]
    }
    }
  else if(all_significant_neg_GIs[i,'query_type']=='min')
  {
    all_significant_neg_GIs[i,'SMF_library'] = supp_table_2$`HAP1_Mutant fitness_minimal (log2 fold change)`[which(supp_table_2$`Gene Symbol`== all_significant_neg_GIs[i,"library_gene"])]
    if(all_significant_neg_GIs[i,"query_gene"] %in% supp_table_2$`Gene Symbol`)
    {
      all_significant_neg_GIs[i,'SMF_query'] = supp_table_2$`HAP1_Mutant fitness_minimal (log2 fold change)`[which(supp_table_2$`Gene Symbol`== all_significant_neg_GIs[i,"query_gene"])]
    }
    }
}

# for each negative GI, calculate metric(qGI_score + single mutant fitness of library gene)
all_significant_neg_GIs$`metric(qGI_score+SMF_library)`= all_significant_neg_GIs$qGI_score + all_significant_neg_GIs$SMF_library

write.csv(all_significant_neg_GIs,"./output/SL_information_NE_OMIM.csv",row.names=F)

# generate plots no. of unique genes and gene pairs at different Double mutant fitness score

limit_list = c(0,-1,-1.5)
step_list = c(0.1,0.1,0.05)
output_folder = "./output/"
for(index in c(1,2,3)){
	#t_list = seq(-3, 0, by = 0.1)
	#t_list = seq(-3, -1, by = 0.1)
	#t_list = seq(-3, -1.5, by = 0.05)
	
	t_list = seq(-3, limit_list[index], by = step_list[index])
	gene_count = c()
	pair_count = c()
	for(t in t_list){
	  gene_count = c(gene_count,length(unique(all_significant_neg_GIs$library_gene[which(all_significant_neg_GIs$`metric(qGI_score+SMF_library)`< (t))])))
	  pair_count = c(pair_count,length(all_significant_neg_GIs$library_gene[which(all_significant_neg_GIs$`metric(qGI_score+SMF_library)`< (t))]))
	}
	all_data = data.frame('threshold' = t_list, 'gene count' = gene_count, 'pair.count' = pair_count)

	plot1 <- ggplot() +
	  geom_line(data = all_data, aes(x = as.numeric(threshold), y=gene.count),colour = 'grey',size=2)+
	  xlab("") +
	  ylab("unique genes") +
	  theme_bw() +
	  theme(
		axis.text = element_text(family="ArialMT",size = 30,colour = "black"),
		legend.title = element_text(family="ArialMT",size = 25,colour = "black"),
		legend.text = element_text(family="ArialMT",size = 25,colour = "black"),
		axis.title = element_text(family="ArialMT",size = 25,colour = "black"),
		panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size = 1)
	  )

	plot2 <- ggplot() +
	  geom_line(data = all_data, aes(x = as.numeric(threshold), y=pair.count),colour = 'grey',size=2)+
	  xlab("Double mutant fitness") +
	  ylab("gene pairs") +
	  theme_bw() +
	  theme(
		axis.text = element_text(family="ArialMT",size = 30,colour = "black"),
		legend.title = element_text(family="ArialMT",size = 25,colour = "black"),
		legend.text = element_text(family="ArialMT",size = 25,colour = "black"),
		axis.title = element_text(family="ArialMT",size = 25,colour = "black"),
		panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size = 1)
	  )

	
	ggsave(file.path(output_folder, paste('NE_OMIM_SL_3_',abs(limit_list[index]),'.pdf',sep='')), arrangeGrob(plot1, plot2),height = 15 , width = 10)
}