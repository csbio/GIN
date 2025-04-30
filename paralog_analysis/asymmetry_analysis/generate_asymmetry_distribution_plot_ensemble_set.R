
library("tibble", lib.loc="../../local/R/library/")
library("readxl", lib.loc="../../local/R/library/")
library("dplyr", lib.loc="../../local/R/library/")
library("reshape2", lib.loc="../../local/R/library/")
library("ggplot2", lib.loc="../../local/R/library/")
library("ggthemes", lib.loc="../../local/R/library/")
library("gridExtra", lib.loc="../../local/R/library/")

output_folder = "./"

#load paralog gene pairs list
ohnolog = readRDS(file='./ensembl_ohnolog_pairs_complete.rds')
ohnolog = ohnolog[!duplicated(ohnolog$sorted_pair),] # keep only one occurence of each pair

#filter based on expression
supp_table_2_file = "./Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx"
supp_table_2 = read_excel(supp_table_2_file)
supp_table_sub = supp_table_2[which(supp_table_2$`HAP1_mRNA expression [log2(TPM+1)]`>1),]

ohnolog_sub_both_exp = ohnolog[which(ohnolog$name1 %in% supp_table_sub$`Gene Symbol` & ohnolog$name2 %in% supp_table_sub$`Gene Symbol`),]

#filter based on sequence id
ohnolog_sub_both_exp = ohnolog_sub_both_exp[which(ohnolog_sub_both_exp$ident_max > 50),]

#degree information non-redundant library
degree_data = read.csv('./qGI_degree_auc_ess.txt',sep = '\t')

ohnolog_sub_both_exp = ohnolog_sub_both_exp[which(ohnolog_sub_both_exp$name1 %in% degree_data$gene & ohnolog_sub_both_exp$name2 %in% degree_data$gene),]

ohnolog_sub_both_exp_ = ohnolog_sub_both_exp[,c('name1','name2','sorted_pair')]

# get negative interaction degree of each gene pair 
ohnolog_sub_both_exp_$qGI_degree_std_gene1 = degree_data$qGI_degree_std_neg[match(ohnolog_sub_both_exp_$name1,degree_data$gene)]
ohnolog_sub_both_exp_$qGI_degree_std_gene2 = degree_data$qGI_degree_std_neg[match(ohnolog_sub_both_exp_$name2,degree_data$gene)]
# calculate total interaction degree of each gene pair
ohnolog_sub_both_exp_$total_degree = ohnolog_sub_both_exp_$qGI_degree_std_gene1 + ohnolog_sub_both_exp_$qGI_degree_std_gene2
# set the lower of the two degrees to numerator
numerator = if_else(ohnolog_sub_both_exp_$qGI_degree_std_gene1 < ohnolog_sub_both_exp_$qGI_degree_std_gene2, ohnolog_sub_both_exp_$qGI_degree_std_gene2, ohnolog_sub_both_exp_$qGI_degree_std_gene1)
# set the higher of the two degrees to denominator
denominator = if_else(ohnolog_sub_both_exp_$qGI_degree_std_gene1 < ohnolog_sub_both_exp_$qGI_degree_std_gene2, ohnolog_sub_both_exp_$qGI_degree_std_gene1, ohnolog_sub_both_exp_$qGI_degree_std_gene2)
# calculate random degree ratio
ohnolog_sub_both_exp_$delta_gene1_gene2 = numerator/denominator

#Filter out pairs with total degree < 4
ohnolog_sub_both_exp_ = ohnolog_sub_both_exp_[which(ohnolog_sub_both_exp_$total_degree>4),]

highest_cutoff = 30
# ratio INF : set to 30 (the bin with ratio>30)
ohnolog_sub_both_exp_$delta_gene1_gene2[!is.finite(ohnolog_sub_both_exp_$delta_gene1_gene2)] <- highest_cutoff
# reduce the right tail - clamp all >30 values to 30
ohnolog_sub_both_exp_$delta_gene1_gene2[ohnolog_sub_both_exp_$delta_gene1_gene2>highest_cutoff] <- highest_cutoff

delta_df = data.frame(ohnolog_sub_both_exp_$delta_gene1_gene2)
df_freq = data.frame(table(delta_df))

# gene pairs with ratio of more than 1/15
sum(df_freq$Freq[which(as.numeric(df_freq$delta_df)>15)])

####generate background null model ###################
#background probability
set.seed(10)
random_delta = c()
random_degree_1_list = c()
random_degree_2_list = c()

for(index in c(1:nrow(ohnolog_sub_both_exp_)))
{
  total_degree = ohnolog_sub_both_exp_$qGI_degree_std_gene1[index] + ohnolog_sub_both_exp_$qGI_degree_std_gene2[index] # get actual total degree of a paralog pair
  random_degree_1 = rbinom(n = 1,size = total_degree,prob = .5) # generate a random degree (representing degree of one gene) with 50% probability from a binomial distribution with upper limit of total degree
  random_degree_2 = total_degree - random_degree_1 # calculate random degree of the second gene
  
  random_degree_1_list = c(random_degree_1_list,random_degree_1)
  random_degree_2_list = c(random_degree_2_list,random_degree_2)
  
  numerator = if_else(random_degree_1 < random_degree_2, random_degree_2, random_degree_1) # set the lower of the two random degrees to numerator
  denominator = if_else(random_degree_1 < random_degree_2, random_degree_1, random_degree_2) # set the higher of the two random degrees to denominator
  random_delta = c(random_delta,numerator/denominator) # calculate random degree ratio and append to list
  
}
background_df = data.frame(random_degree_1_list,random_degree_2_list,random_delta) 
# calculate total random degree of each gene pair
background_df$total_degree = background_df$random_degree_1_list + background_df$random_degree_2_list
#Filter out pairs with total degree < 4
background_df = background_df[which(background_df$total_degree>4),]

# ratio INF : set to 30 (the bin with ratio>30)
background_df$random_delta[!is.finite(background_df$random_delta)] <- highest_cutoff

# reduce the right tail - clamp all >30 values to 30
background_df$random_delta[background_df$random_delta>highest_cutoff] <- highest_cutoff

###generate distribution plot#################

plot1 <- ggplot() +
  geom_density(data = background_df, aes(x = random_delta, y=after_stat(density)), fill=rgb(210, 210, 210,maxColorValue = 255), color = NA,alpha = .5,kernel = "gaussian",adjust = 5) +
  geom_histogram(data = delta_df, aes(x = ohnolog_sub_both_exp_.delta_gene1_gene2, y=after_stat(density),fill=after_stat(x)), bins = 20,alpha = .5) +
  xlab("ratio") +
  ylab("Density") +
  scale_x_continuous(breaks=seq(0,30,5),labels=c(0 , '1:5', '1:10', '1:15', '1:20', '1:25', '1:>=30')) +
  
  theme_bw() +
  theme(
    axis.text = element_text(family="ArialMT",size = 30,colour = "black"),
    legend.title = element_text(family="ArialMT",size = 25,colour = "black"),
    legend.text = element_text(family="ArialMT",size = 25,colour = "black"),
    axis.title = element_text(family="ArialMT",size = 25,colour = "black"),
    panel.border = element_blank(), panel.grid.major = element_blank(), legend.position = "none",
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size = 1)
  )

ggsave(file.path(output_folder, paste0("density_negative_GI_ratio_seqid_50_.pdf")), device = cairo_pdf,height = 10 , width = 10 )

