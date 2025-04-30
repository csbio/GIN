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

#####################################################################################
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
#write.csv(overlap_index_info,"./output_library_3_or_query_1_complex/removed_complexes_subset_5_.1_between_version_2.csv",row.names = F)

# list the subset complexes (that has a larger superset complex with overlap index of .3) to remove
#not_keep_1 = unique(overlap_index_info_1$com_2)
#not_keep_2 = unique(overlap_index_info_2$com_2)
not_keep_2 = unique(overlap_index_info$com_2)

# remove the subset complexes					 
#data_qGI_between <- data_qGI_sub[which(!(data_qGI_sub$Row_Complex %in% not_keep_1 & data_qGI_sub$Column_Complex %in% not_keep_2)),]
data_qGI_sub <- data_qGI_sub[which(!(data_qGI_sub$Row_Complex %in% not_keep_2 | data_qGI_sub$Column_Complex %in% not_keep_2)),]

#keep rows with total enrichment<.1
data_qGI_sub <- data_qGI_sub[which(data_qGI_sub$Enrichment_total<.1),]

# save output
write.csv(data_qGI_sub,file.path("./output_library_3_or_query_1_complex/NO_SUBSET_.3_WB_library_3_or_query_1_complex_subset_5_.1_between_version_2.csv"),row.names = F)

#get 
# list of # of positive GIs per module, 
# list of # of negative GIs per module, 
# purity scores (ranges from -1 to 1, 1 being all positive interaction and -1 being all negative interaction), 
#	purity = actual pos. GI in a module / actual all GI in a module
# list of # of total GIs
observed_pos = data_qGI_sub$No_pos_interactions
observed_neg = data_qGI_sub$No_neg_interactions
observed_purity = data_qGI_sub$Purity_score
observed_total = data_qGI_sub$No_total_interactions
prob = sum(observed_pos)/(sum(observed_pos)+sum(observed_neg))

# for each module, 
# generate a random positive GI given the # of total GIs and the background probability as weight, from a binomial distribution 
# and calculate the random (background) purity score.
set.seed(10)
random_pos = c()
random_purity = c()
for(index in c(1:nrow(data_qGI_sub)))
{
  pos_ = rbinom(n = 1,size = observed_total[index],prob = prob)
  random_pos = c(random_pos,pos_)
  random_purity = c(random_purity,pos_/observed_total[index])
}

# normalize random (background) purity scores to range [-1,1]
background = rescale(random_purity, to = c(-1, 1))

background_df = data.frame(background)
purity_df = data.frame(observed_purity)

#generate distribution barplot of purity scores against random (background) purity density plot
plot1 <- ggplot() +
  geom_density(data = background_df, aes(x = background, y=after_stat(density)), fill=rgb(210, 210, 210,maxColorValue = 255), color = NA,alpha = .5,kernel = "gaussian",adjust = 1.5) +
  geom_histogram(data = purity_df, aes(x = observed_purity, y=after_stat(density),fill=after_stat(x)), bins = 20) +
  xlab("Purity") +
  #xlim(-1, 1) +
  #ylim(0,1) +
  ylab("Density") +
  #scale_fill_gradient2(low=c(30, 70, 145),  mid = c(200, 200, 200), high=c(250, 220, 0))+
  scale_fill_gradient2(low=rgb(30, 70, 145,maxColorValue = 255),  mid = rgb(200, 200, 200,maxColorValue = 255), high=rgb(250, 220, 0, maxColorValue = 255))+
  
  #labs(color='') +
  theme_bw() +
  theme(
    #axis.title.x=element_blank(),axis.title.y=element_blank(), axis.text.y = element_blank(),
    axis.text = element_text(family="ArialMT",size = 30,colour = "black"),
    legend.title = element_text(family="ArialMT",size = 25,colour = "black"),
    legend.text = element_text(family="ArialMT",size = 25,colour = "black"),
    axis.title = element_text(family="ArialMT",size = 25,colour = "black"),
    panel.border = element_blank(), panel.grid.major = element_blank(), legend.position = "none",
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size = 1)
  )

# output file
ggsave(file.path('./output_library_3_or_query_1_complex/NO_SUBSET_.3_WB_library_3_or_query_1_complex_subset_5_.1_between_version_2.pdf'),plot =plot1,height = 14 , width = 14)

#####################################################################################
#load within complex analysis file: WB_library_3_or_query_1_complex_subset_5_.1_within_.csv 
data_qGI <- read.csv('./output_library_3_or_query_1_complex/WB_library_3_or_query_1_complex_subset_5_within.csv',
                     header = TRUE, stringsAsFactors = FALSE)
#remove same name complex
data_qGI_sub = data_qGI %>% group_by(Row_Complex) %>% top_n(1, No_of_interaction_tested) %>% top_n(1, No_total_interactions)
interesting_complex_set = unique(data_qGI_sub$Row_Complex)
					 
# get gene overlap index information between complexes from subset_5_.1_within_overlap_index_2.csv
overlap_index_info = read.csv('./output_library_3_or_query_1_complex/subset_2_within_overlap_index_version_2.csv',stringsAsFactors = F)
#Get information on complexes for which overlap_index is greater than or equal to .3
# Overlap Index: (size of intersection of two complexes) / (minimum of the two complex sizes)
overlap_index_info = overlap_index_info[which(overlap_index_info$overlap_index>=.3),]
overlap_index_info = overlap_index_info[which((overlap_index_info$com_1 %in% interesting_complex_set) & (overlap_index_info$com_2 %in% interesting_complex_set) ),]

# list the subset complexes (that has a larger superset complex with overlap index of .3) to remove
not_keep = unique(overlap_index_info$com_2)


# remove the subset complexes
data_qGI_sub <- data_qGI_sub[which(!(data_qGI_sub$Row_Complex %in% not_keep)),]

#keep rows with total enrichment<.1
data_qGI_sub <- data_qGI_sub[which(data_qGI_sub$Enrichment_total<.1),]

#save output 
write.csv(data_qGI_sub,file.path("./output_library_3_or_query_1_complex/NO_SUBSET_.3_WB_library_3_or_query_1_complex_subset_5_.1_within_version_2.csv"),row.names = F)

#get 
# list of # of positive GIs per module, 
# list of # of negative GIs per module, 
# purity scores (ranges from -1 to 1, 1 being all positive interaction and -1 being all negative interaction), 
#	purity = actual pos. GI in a module / actual all GI in a module
# list of # of total GIs
observed_pos = data_qGI_sub$No_pos_interactions
observed_neg = data_qGI_sub$No_neg_interactions
observed_purity = data_qGI_sub$Purity_score
observed_total = data_qGI_sub$No_total_interactions
prob = sum(observed_pos)/(sum(observed_pos)+sum(observed_neg))

# for each module, 
# generate a random positive GI given the # of total GIs and the background probability as weight, from a binomial distribution 
# and calculate the random (background) purity score.
set.seed(10)
random_pos = c()
random_purity = c()
for(index in c(1:nrow(data_qGI_sub)))
{
  pos_ = rbinom(n = 1,size = observed_total[index],prob = prob)
  random_pos = c(random_pos,pos_)
  random_purity = c(random_purity,pos_/observed_total[index])
}

# normalize random (background) purity scores to range [-1,1]
background = rescale(random_purity, to = c(-1, 1))

background_df = data.frame(background)
purity_df = data.frame(observed_purity)

#generate distribution barplot of purity scores against random (background) purity density plot
plot1 <- ggplot() +
  geom_density(data = background_df, aes(x = background, y=after_stat(density)), fill=rgb(210, 210, 210,maxColorValue = 255), color = NA,alpha = .5,kernel = "gaussian",adjust = 1.5) +
  geom_histogram(data = purity_df, aes(x = observed_purity, y=after_stat(density),fill=after_stat(x)), bins = 20) +
  xlab("Purity") +
  #xlim(-1, 1) +
  #ylim(0,1) +
  ylab("Density") +
  #scale_fill_gradient2(low=c(30, 70, 145),  mid = c(200, 200, 200), high=c(250, 220, 0))+
  scale_fill_gradient2(low=rgb(30, 70, 145,maxColorValue = 255),  mid = rgb(200, 200, 200,maxColorValue = 255), high=rgb(250, 220, 0, maxColorValue = 255))+
  
  #labs(color='') +
  theme_bw() +
  theme(
    #axis.title.x=element_blank(),axis.title.y=element_blank(), axis.text.y = element_blank(),
    axis.text = element_text(family="ArialMT",size = 30,colour = "black"),
    legend.title = element_text(family="ArialMT",size = 25,colour = "black"),
    legend.text = element_text(family="ArialMT",size = 25,colour = "black"),
    axis.title = element_text(family="ArialMT",size = 25,colour = "black"),
    panel.border = element_blank(), panel.grid.major = element_blank(), legend.position = "none",
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size = 1)
  )

# output file

ggsave(file.path('./output_library_3_or_query_1_complex/NO_SUBSET_.3_WB_library_3_or_query_1_complex_subset_5_.1_within_version_2.pdf'),plot =plot1,height = 14 , width = 14)


