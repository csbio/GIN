library("tibble", lib.loc="../../local/R/library/")
library("readxl", lib.loc="../../local/R/library/")
library("dplyr", lib.loc="../../local/R/library/")

###paralog##################################################################################
#load paralog gene pairs list
ohnolog = readRDS(file='./ensembl_ohnolog_pairs_complete.rds')
ohnolog = ohnolog[!duplicated(ohnolog$sorted_pair),]

# filter gene pairs based on expression
supp_table_2_file = "./Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx"
supp_table_2 = read_excel(supp_table_2_file)
supp_table_sub = supp_table_2[which(supp_table_2$`HAP1_mRNA expression [log2(TPM+1)]`>1),]

ohnolog_sub_both_exp = ohnolog[which(ohnolog$name1 %in% supp_table_sub$`Gene Symbol` & ohnolog$name2 %in% supp_table_sub$`Gene Symbol`),]

# filter gene pairs based on sequence id
ohnolog_sub_both_exp = ohnolog_sub_both_exp[which(ohnolog_sub_both_exp$ident_max > 50),]

# get degree information of non-redundant library genes
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

#Keep gene pairs with total degree > 4
ohnolog_sub_both_exp_ = ohnolog_sub_both_exp_[which(ohnolog_sub_both_exp_$total_degree>4),]

# set degree equal to a cutoff if above a cut-off 
highest_cutoff = 30
ohnolog_sub_both_exp_$delta_gene1_gene2[!is.finite(ohnolog_sub_both_exp_$delta_gene1_gene2)] <- highest_cutoff
ohnolog_sub_both_exp_$delta_gene1_gene2[ohnolog_sub_both_exp_$delta_gene1_gene2>highest_cutoff] <- highest_cutoff

delta_df = data.frame(ohnolog_sub_both_exp_$delta_gene1_gene2)

average_real_delta = mean(delta_df$ohnolog_sub_both_exp_.delta_gene1_gene2)

# generate random background degree ratio

average_null_delta = c()
for(i in 1:1000) # iterate for 1000 times
{
print(i)
  random_delta = c()
	random_degree_1_list = c()
	random_degree_2_list = c()

	for(index in c(1:nrow(ohnolog_sub_both_exp_))) #generate random degrees for all paralog gene pairs
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
  background_df$total_degree = background_df$random_degree_1_list + background_df$random_degree_2_list # calculate total random degree of each gene pair
  background_df = background_df[which(background_df$total_degree>4),] # keep gene pairs with total degree > 4
  # set random degree equal to a cutoff if above a cut-off 
  background_df$random_delta[!is.finite(background_df$random_delta)] <- highest_cutoff
  background_df$random_delta[background_df$random_delta>highest_cutoff] <- highest_cutoff 

  average_null_delta = c(average_null_delta,mean(background_df$random_delta)) # calculate average random ratio across gene pairs and append to list
}

write.csv(average_null_delta,"./average_null_negative_GI_ratio_seqid_50.csv",row.names=F)

print(average_real_delta) 
# https://rdrr.io/bioc/qvalue/man/empPvals.html
# https://rdrr.io/bioc/qvalue/src/R/empPvals.R 

null_delta = read.csv('./average_null_negative_GI_ratio_seqid_50.csv')
print(mean(null_delta$x)) 
real_delta = c(average_real_delta) 
result = empPvals(real_delta, null_delta$x, pool = TRUE)

print(result)

p_fun <- ecdf(null_delta$x)
emp_p <- 1-unlist(lapply(real_delta, p_fun))
print(emp_p)