options(java.parameters = "- Xmx1024m")
library("withr", lib.loc="../../local/R/library/")
library("farver", lib.loc="../../local/R/library/")
library("labeling", lib.loc="../../local/R/library/")
library("digest", lib.loc="../../local/R/library/")
library("readxl", lib.loc="../../local/R/library/")
library("dplyr", lib.loc="../../local/R/library/")
library("reshape2", lib.loc="../../local/R/library/")

data_qGI_nRed <- read.table('./qGI_20211111_nRed.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE) 
data_fdr_nRed <- read.table('./gi_FDR_20211111_nRed.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE) 

data_qGI_nRed_melt <- setNames(melt(as.matrix(data_qGI_nRed)), c('library_gene', ' query_name', 'qGI_score'))
data_fdr_nRed_melt <- setNames(melt(as.matrix(data_fdr_nRed)), c('library_gene', ' query_name', 'FDR'))

all_GIs = merge(data_qGI_nRed_melt,data_fdr_nRed_melt,by = c('library_gene', ' query_name'))

split_screen_name<-function(x)
{
  x1<-strsplit(x,"_")[[1]][1] 
  return( x1)
}
all_GIs$query_gene <- unlist(lapply(as.character(all_GIs$` query_name`),split_screen_name) )
all_GIs <- all_GIs %>% mutate_if(is.factor, as.character)
all_GIs = all_GIs[-which(all_GIs$` query_name` =='TAPT1_375_rich'),]

all_GIs = all_GIs %>%
  rowwise() %>%
  mutate(sorted_pair = paste(unlist(sort(c(as.character(library_gene),query_gene))), collapse = "_"))
  
write.csv(all_GIs$sorted_pair,"./all_GIs_gene_pairs.csv",row.names=F)

all_significant_positive_GIs = all_GIs[which(all_GIs$FDR<.1 & all_GIs$qGI_score>.3),]

write.csv(all_significant_positive_GIs$sorted_pair,"./all_significant_positive_GIs_gene_pairs.csv",row.names=F)


all_significant_negative_GIs = all_GIs[which(all_GIs$FDR<.1 & all_GIs$qGI_score<(-.3)),]

write.csv(all_significant_negative_GIs$sorted_pair,"./all_significant_negative_GIs_gene_pairs.csv",row.names=F)
