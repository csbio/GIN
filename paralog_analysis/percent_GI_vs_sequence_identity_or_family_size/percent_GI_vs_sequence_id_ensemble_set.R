library("withr", lib.loc="../../local/R/library/")
library("farver", lib.loc="../../local/R/library/")
library("labeling", lib.loc="../../local/R/library/")
library("digest", lib.loc="../../local/R/library/")
library("readxl", lib.loc="../../local/R/library/")
library("dplyr", lib.loc="../../local/R/library/")
library("reshape2", lib.loc="../../local/R/library/")
library("ggplot2", lib.loc="../../local/R/library/")
library("ggthemes", lib.loc="../../local/R/library/")
library("gridExtra", lib.loc="../../local/R/library/")


output_folder = "./"

# load significant positive and negative GI gene pairs
all_significant_positive_GIs = read.csv("./all_significant_positive_GIs_gene_pairs.csv",header=T)
all_significant_negative_GIs = read.csv("./all_significant_negative_GIs_gene_pairs.csv",header=T)
all_GIs = read.csv("./all_GIs_gene_pairs.csv",header=T)

# calculate background percent of positive and negative GI gene pairs
background_pos = nrow(all_significant_positive_GIs)/nrow(all_GIs)*100 
background_neg = nrow(all_significant_negative_GIs)/nrow(all_GIs)*100

# load ensemble paralog information
ohnolog = readRDS(file='./ensembl_ohnolog_pairs_complete.rds')
ohnolog = ohnolog[!duplicated(ohnolog$sorted_pair),] # remove duplicate pairs

# filter screened gene pair list based on expression - keep gene pairs with only log2(TPM+1) expression > 1
supp_table_2_file = "./Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx"
supp_table_2 = read_excel(supp_table_2_file)
supp_table_sub = supp_table_2[which(supp_table_2$`HAP1_mRNA expression [log2(TPM+1)]`>1),]

# get paralogs where neither gene meets the expression cut-off
ohnolog_sub_no_exp = ohnolog[- (which(ohnolog$name1 %in% supp_table_sub$`Gene Symbol` | ohnolog$name2 %in% supp_table_sub$`Gene Symbol`)),]
# get paralogs where both gene meets the expression cut-off
ohnolog_sub_both_exp = ohnolog[which(ohnolog$name1 %in% supp_table_sub$`Gene Symbol` & ohnolog$name2 %in% supp_table_sub$`Gene Symbol`),]
# get paralogs where only one gene meets the expression cut-off
ohnolog_sub_one_exp = ohnolog[which(ohnolog$name1 %in% supp_table_sub$`Gene Symbol` | ohnolog$name2 %in% supp_table_sub$`Gene Symbol`),]
ohnolog_sub_one_exp = ohnolog_sub_one_exp[which(ohnolog_sub_one_exp$sorted_pair %in% setdiff(ohnolog_sub_one_exp$sorted_pair,ohnolog_sub_both_exp$sorted_pair)),]

expression_filter = 'no_exp_vs_both_exp_vs_one_exp_' #'no_exp' 'both' 'atleast_1'

#############################################################
id_list_left = seq(20, 90, by=10)

exp_df = data.frame(matrix(nrow=0,ncol=6))

# iterate for sequence id bins 
for(id_index in 1:length(id_list_left))
{
	#for each category of paralogs, keep gene pairs that were screened and calculate percentage of positive and negative GIs
	
  ohnolog_all_no_exp = ohnolog_sub_no_exp$sorted_pair[which(ohnolog_sub_no_exp$ident_max>id_list_left[id_index] )]
  ohnolog_screened_no_exp = unique(intersect(ohnolog_all_no_exp, all_GIs$x))
  ohnolog_screened_pos_no_exp = unique(intersect(ohnolog_all_no_exp, all_significant_positive_GIs$x))
  ohnolog_screened_neg_no_exp = unique(intersect(ohnolog_all_no_exp, all_significant_negative_GIs$x))
  percent_pos_no_exp = length(ohnolog_screened_pos_no_exp)/length(ohnolog_screened_no_exp) * 100
  percent_neg_no_exp = length(ohnolog_screened_neg_no_exp)/length(ohnolog_screened_no_exp) * 100
  
  ohnolog_all_both_exp = ohnolog_sub_both_exp$sorted_pair[which(ohnolog_sub_both_exp$ident_max>id_list_left[id_index] )]
  ohnolog_screened_both_exp = unique(intersect(ohnolog_all_both_exp, all_GIs$x))
  ohnolog_screened_pos_both_exp = unique(intersect(ohnolog_all_both_exp, all_significant_positive_GIs$x))
  ohnolog_screened_neg_both_exp = unique(intersect(ohnolog_all_both_exp, all_significant_negative_GIs$x))
  percent_pos_both_exp = length(ohnolog_screened_pos_both_exp)/length(ohnolog_screened_both_exp) * 100
  percent_neg_both_exp = length(ohnolog_screened_neg_both_exp)/length(ohnolog_screened_both_exp) * 100
  
  ohnolog_all_one_exp = ohnolog_sub_one_exp$sorted_pair[which(ohnolog_sub_one_exp$ident_max>id_list_left[id_index] )]
  ohnolog_screened_one_exp = unique(intersect(ohnolog_all_one_exp, all_GIs$x))
  ohnolog_screened_pos_one_exp = unique(intersect(ohnolog_all_one_exp, all_significant_positive_GIs$x))
  ohnolog_screened_neg_one_exp = unique(intersect(ohnolog_all_one_exp, all_significant_negative_GIs$x))
  percent_pos_one_exp = length(ohnolog_screened_pos_one_exp)/length(ohnolog_screened_one_exp) * 100
  percent_neg_one_exp = length(ohnolog_screened_neg_one_exp)/length(ohnolog_screened_one_exp) * 100
  
  exp_df = rbind(exp_df,
                 c(
                   paste('>',id_list_left[id_index],
                         ' (',length(ohnolog_screened_no_exp),
                         ',',length(ohnolog_screened_both_exp),
                         ',',length(ohnolog_screened_one_exp),' )',sep = ''),
                  as.numeric(percent_pos_no_exp),
                  as.numeric(percent_neg_no_exp),
                  as.numeric(percent_pos_both_exp),
                  as.numeric(percent_neg_both_exp),
                  as.numeric(percent_pos_one_exp),
                  as.numeric(percent_neg_one_exp)
                  ),stringsAsFactors = FALSE
  )

}

#prepare data for plotting
colnames(exp_df) = c('id_bin',
                     'percent_pos_no_exp', 'percent_neg_no_exp',
                     'percent_pos_both_exp', 'percent_neg_both_exp',
                     'percent_pos_one_exp', 'percent_neg_one_exp'
                     )

exp_df$percent_pos_no_exp[is.nan(as.numeric(exp_df$percent_pos_no_exp))] = 0
exp_df$percent_neg_no_exp[is.nan(as.numeric(exp_df$percent_neg_no_exp))] = 0
exp_df_melt <- melt(exp_df ,  id.vars = 'id_bin', variable.name = 'series')

cols <- c("percent_pos_no_exp" = "orange",
          "percent_neg_no_exp" = "lightblue", 
          "percent_pos_both_exp" = "yellow",
          "percent_neg_both_exp" = "blue",
          "percent_pos_one_exp" = "red",
          "percent_neg_one_exp" = "darkblue"
          )
# generate line plot - percent of positive and negative GI vs sequence ID bin
PVEplot <- ggplot(data = exp_df, aes(x=as.factor(id_bin))) +
  geom_point( aes(y=as.numeric(percent_pos_no_exp)),color='orange',size=4) +
  geom_point(aes(y=as.numeric(percent_neg_no_exp)),color='lightblue',size=4) +
  geom_point(aes(y=as.numeric(percent_pos_both_exp)),color='yellow',size=4) +
  geom_point(aes(y=as.numeric(percent_neg_both_exp)),color='blue',size=4) +
  geom_point(aes(y=as.numeric(percent_pos_one_exp)),color='red',size=4) +
  geom_point(aes(y=as.numeric(percent_neg_one_exp)),color='darkblue',size=4) +

  geom_hline(yintercept = background_pos, color='yellow',linetype=2) +
  geom_hline(yintercept = background_neg, color='blue',linetype=2) +
  
  geom_line(data=exp_df_melt,aes(as.factor(id_bin), as.numeric(value), group=series,colour = factor(series)), )+
  scale_color_manual(values = cols,
                     labels=c('no exp. pos.','no exp. neg.','both exp. pos.','both exp. neg.','one exp. pos.','one exp. neg.')) +

  labs(color = '')+
  xlab("ID bin") +
  ylab("% of GI") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), #axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        legend.title = element_text(family="Sans",size = 20,colour = "black"),
        legend.text = element_text(family="Sans",size = 20,colour = "black"),
        axis.title = element_text(family="Sans",size = 20,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(family="Sans",size = 20,colour = "black"),
        axis.text.x = element_text(family="Sans",size = 20,colour = "black",angle = 90,hjust=0.95,vjust=0.5),
        plot.margin = unit(c(1,1,1,1), "cm"),
        axis.ticks.length=unit(0.1,"inch")
  )

ggsave(file.path(output_folder, paste0(expression_filter,"_exp_id_open_bin_dot_plot_.pdf")), device = cairo_pdf,height = 10 , width = 15 )
