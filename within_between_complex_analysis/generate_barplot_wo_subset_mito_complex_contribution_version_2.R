library("withr", lib.loc="../../local/R/library/")
library("farver", lib.loc="../../local/R/library/")
library("labeling", lib.loc="../../local/R/library/")
library("digest", lib.loc="../../local/R/library/")
library("ggplot2", lib.loc="../../local/R/library/")
library("ggthemes", lib.loc="../../local/R/library/")
library("RColorBrewer", lib.loc="../../local/R/library/")
library("gridExtra", lib.loc="../../local/R/library/")
library("scales", lib.loc="../../local/R/library/")
library("dplyr", lib.loc="../../local/R/library/")

######################################################################################

# get complexes in which the fraction of mitochondrial genes is more than .5
# mitochondrial_fraction_complexes.csv contains information regarding fraction of mitochondrial genes in the complexes
mito_info = read.csv('./output_library_3_or_query_1_complex/mitochondrial_fraction_complexes.csv',stringsAsFactors = F)
mito_complex <- mito_info$complex[which(mito_info$fraction_mito>0.5)]

###between###################################################################################

#load between complex pair analysis file: WB_library_3_or_query_1_complex_subset_5_.1_between_.csv 
data_qGI <- read.csv('./output_library_3_or_query_1_complex/WB_library_3_or_query_1_complex_subset_5_between.csv',
                     header = TRUE, stringsAsFactors = FALSE)
					 
#remove same name complex
data_qGI_sub = data_qGI %>% group_by(Row_Complex,Column_Complex) %>% top_n(1, No_of_interaction_tested) %>% top_n(1, No_total_interactions)

#interesting_complex_set_1 = unique(c(data_qGI_sub$Row_Complex))
#interesting_complex_set_2 = unique(c(data_qGI_sub$Column_Complex))
interesting_complex_set = unique(c(data_qGI_sub$Row_Complex,data_qGI_sub$Column_Complex))
					
# get gene overlap index information between complexes from subset_5_.1_between_overlap_index_2.csv
overlap_index_info = read.csv('./output_library_3_or_query_1_complex/subset_5_between_overlap_index_version_2.csv',stringsAsFactors = F)
#Get information on complexes for which overlap_index is greater than or equal to .3
# Overlap Index: (size of intersection of two complexes) / (minimum of the two complex sizes)
overlap_index_info = overlap_index_info[which(overlap_index_info$overlap_index>=.3),]
#overlap_index_info_1 = overlap_index_info[which((overlap_index_info$com_1 %in% interesting_complex_set_1) & (overlap_index_info$com_2 %in% interesting_complex_set_1)),]
#overlap_index_info_2 = overlap_index_info[which((overlap_index_info$com_1 %in% interesting_complex_set_2) & (overlap_index_info$com_2 %in% interesting_complex_set_2)),]
overlap_index_info = overlap_index_info[which((overlap_index_info$com_1 %in% interesting_complex_set) & (overlap_index_info$com_2 %in% interesting_complex_set) ),]

# list the subset complexes (that has a larger superset complex with overlap index of .3) to remove
#not_keep_1 = unique(overlap_index_info_1$com_2)
#not_keep_2 = unique(overlap_index_info_2$com_2)
not_keep_2 = unique(overlap_index_info$com_2)

# remove the subset complexes					 
#data_qGI_between <- data_qGI_sub[which(!(data_qGI_sub$Row_Complex %in% not_keep_1 & data_qGI_sub$Column_Complex %in% not_keep_2)),]
data_qGI_between <- data_qGI_sub[which(!(data_qGI_sub$Row_Complex %in% not_keep_2 | data_qGI_sub$Column_Complex %in% not_keep_2)),]

# remove mitochondrial complexes
data_qGI_between <- data_qGI_between[which(!(data_qGI_between$Row_Complex %in% mito_complex | data_qGI_between$Column_Complex %in% mito_complex)),]

###within####################################################################

#load within complex analysis file: WB_library_3_or_query_1_complex_subset_5_.1_within_.csv 
data_qGI <- read.csv('./output_library_3_or_query_1_complex/WB_library_3_or_query_1_complex_subset_5_within.csv',
                     header = TRUE, stringsAsFactors = FALSE)

#remove same name complex
data_qGI_sub = data_qGI %>% group_by(Row_Complex) %>% top_n(1, No_of_interaction_tested) %>% top_n(1, No_total_interactions)
interesting_complex_set = unique(data_qGI_sub$Row_Complex)
					 
# get gene overlap index information between complexes from subset_5_.1_within_overlap_index_2.csv
overlap_index_info = read.csv('./output_library_3_or_query_1_complex/subset_5_within_overlap_index_version_2.csv',stringsAsFactors = F)
#Get information on complexes for which overlap_index is greater than or equal to .3
# Overlap Index: (size of intersection of two complexes) / (minimum of the two complex sizes)
overlap_index_info = overlap_index_info[which(overlap_index_info$overlap_index>=.3),]
overlap_index_info = overlap_index_info[which((overlap_index_info$com_1 %in% interesting_complex_set) & (overlap_index_info$com_2 %in% interesting_complex_set) ),]

# list the subset complexes (that has a larger superset complex with overlap index of .3) to remove
not_keep = unique(overlap_index_info$com_2)

# remove the subset complexes					 
data_qGI_within <- data_qGI_sub[which(!(data_qGI_sub$Row_Complex %in% not_keep)),]

# remove mitochondrial complexes
data_qGI_within <- data_qGI_within[which(!(data_qGI_within$Row_Complex %in% mito_complex | data_qGI_within$Column_Complex %in% mito_complex)),]

## store to file 

write.csv(data_qGI_within,file.path("./output_library_3_or_query_1_complex/NO_MITO_.5_NO_SUBSET_.3_WB_library_3_or_query_1_complex_subset_5_within_version_2.csv"),row.names = F)
write.csv(data_qGI_between,file.path("./output_library_3_or_query_1_complex/NO_MITO_.5_NO_SUBSET_.3_WB_library_3_or_query_1_complex_subset_5_between_version_2.csv"),row.names = F)

#######################################################################

# generate dataframe with GI stats
# percentage of complexes with non-zero positive GIs with pos. enrichment score smaller than .1
# percentage of complexes with non-zero negative GIs with neg. enrichment score smaller than .1
# percentage of complexes with non-zero total GIs with total enrichment score smaller than .1
w_list = data.frame(c('positive','negative','all'))
colnames(w_list) = c('interaction')
w_list$w_complex_percentage = c(sum(data_qGI_within$No_pos_interactions>0 & data_qGI_within$Enrichment_pos<.1)*100/nrow(data_qGI_within),
                                sum(data_qGI_within$No_neg_interactions>0 & data_qGI_within$Enrichment_neg<.1)*100/nrow(data_qGI_within),
                                sum(data_qGI_within$No_total_interactions>0 & data_qGI_within$Enrichment_total<.1)*100/nrow(data_qGI_within)
)


w_list$b_complex_percentage = c(sum(data_qGI_between$No_pos_interactions>0 & data_qGI_between$Enrichment_pos<.1)*100/nrow(data_qGI_between),
                                sum(data_qGI_between$No_neg_interactions>0 & data_qGI_between$Enrichment_neg<.1)*100/nrow(data_qGI_between),
                                sum(data_qGI_between$No_total_interactions>0 & data_qGI_between$Enrichment_total<.1)*100/nrow(data_qGI_between)
)

# generate within complex GI stat. barplot
plot1 <- ggplot(data=w_list) +
  geom_col(aes(x = interaction, y = w_complex_percentage),  alpha = 1) +
  ggtitle("Within complex") +
  xlab("") +
  ylab("Percentage") +
  labs(color='') +
  theme_bw() +
  coord_flip() +
  theme(
    #axis.title.x=element_blank(),axis.title.y=element_blank(), axis.text.y = element_blank(),
    axis.text = element_text(family="ArialMT",size = 30,colour = "black"),
    legend.title = element_text(family="ArialMT",size = 25,colour = "black"),
    legend.text = element_text(family="ArialMT",size = 25,colour = "black"),
    axis.title = element_text(family="ArialMT",size = 25,colour = "black"),
    panel.border = element_blank(), panel.grid.major = element_blank(), legend.position = "none",
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size = 1)
  )

# generate between complex GI stat. barplot
plot2 <- ggplot(data=w_list) +
  geom_col(aes(x = interaction, y = b_complex_percentage),  alpha = 1) +
  ggtitle("Between complex") +
  xlab("") +
  ylab("Percentage") +
  labs(color='') +
  theme_bw() +
  coord_flip() +
  theme(
    #axis.title.x=element_blank(),axis.title.y=element_blank(), axis.text.y = element_blank(),
    axis.text = element_text(family="ArialMT",size = 30,colour = "black"),
    legend.title = element_text(family="ArialMT",size = 25,colour = "black"),
    legend.text = element_text(family="ArialMT",size = 25,colour = "black"),
    axis.title = element_text(family="ArialMT",size = 25,colour = "black"),
    panel.border = element_blank(), panel.grid.major = element_blank(), legend.position = "none",
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size = 1)
  )

plot3 <- grid.arrange(plot1, plot2,nrow = 2)
ggsave(file.path('./output_library_3_or_query_1_complex/NO_MITO_.5_NO_SUBSET_.3_complex_percentage_bar_plot_wb_5_version_2.pdf'),plot =plot3)
