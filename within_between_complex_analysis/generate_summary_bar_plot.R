

library("withr", lib.loc="../../local/R/library/")
library("farver", lib.loc="../../local/R/library/")
library("labeling", lib.loc="../../local/R/library/")
library("digest", lib.loc="../../local/R/library/")
library("ggplot2", lib.loc="../../local/R/library/")
library("ggthemes", lib.loc="../../local/R/library/")
library("RColorBrewer", lib.loc="../../local/R/library/")
library("gridExtra", lib.loc="../../local/R/library/")
library("scales", lib.loc="../../local/R/library/")

################################################################################
#load within complex pair analysis file: WB_library_3_or_query_1_complex_subset_5_within.csv 
data_qGI_within <- read.csv('./output_library_3_or_query_1_complex/WB_library_3_or_query_1_complex_subset_5_within.csv',
                     header = TRUE, stringsAsFactors = FALSE)

#load between complex pair analysis file: WB_library_3_or_query_1_complex_subset_5_between.csv 
data_qGI_between <- read.csv('./output_library_3_or_query_1_complex/WB_library_3_or_query_1_complex_subset_5_between.csv',
                     header = TRUE, stringsAsFactors = FALSE)				

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

ggsave(file.path('./output_library_3_or_query_1_complex/complex_percentage_bar_plot_wb_5_.pdf'),plot =plot3)
