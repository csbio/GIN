
library("withr", lib.loc="../../local/R/library/")
library("farver", lib.loc="../../local/R/library/")
library("labeling", lib.loc="../../local/R/library/")
library("digest", lib.loc="../../local/R/library/")
library("ggplot2", lib.loc="../../local/R/library/")
library("ggthemes", lib.loc="../../local/R/library/")
library("RColorBrewer", lib.loc="../../local/R/library/")
#library("rstan", lib.loc="../../local/R/library/")
library("scales", lib.loc="../../local/R/library/")

################################################################################
#load between complex pair analysis file: WB_library_3_or_query_1_complex_subset_5_.1_between.csv 
data_qGI <- read.csv('./output_library_3_or_query_1_complex/WB_library_3_or_query_1_complex_subset_5_.1_between.csv',
                     header = TRUE, stringsAsFactors = FALSE)

#get 
# list of # of positive GIs per module, 
# list of # of negative GIs per module, 
# purity scores (ranges from -1 to 1, 1 being all positive interaction and -1 being all negative interaction), 
#	purity = actual pos. GI in a module / actual all GI in a module
# list of # of total GIs
observed_pos = data_qGI$No_pos_interactions
observed_neg = data_qGI$No_neg_interactions
observed_purity = data_qGI$Purity_score
observed_total = data_qGI$No_total_interactions

# background probability = sum of all positive GIs across all modules / sum of all GIs across all modules
prob = sum(observed_pos)/(sum(observed_pos)+sum(observed_neg))

set.seed(10)

#x = rbinom(n = nrow(data_qGI),size = 100,prob = prob)
#background = rescale(x, to = c(-1, 1), from = c(0,100))

# for each module, 
# generate a random positive GI given the # of total GIs and the background probability as weight, from a binomial distribution 
# and calculate the random (background) purity score.
random_pos = c()
random_purity = c()
for(index in c(1:nrow(data_qGI)))
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

# output files

ggsave(file.path('./output_library_3_or_query_1_complex/', 'WB_library_3_or_query_1_complex_subset_5_.1_between_final.pdf'),plot =plot1,height = 14 , width = 14)

###################################################################################
#load within complex analysis file: WB_library_3_or_query_1_complex_subset_5_.1_within.csv
data_qGI <- read.csv('./output_library_3_or_query_1_complex/WB_library_3_or_query_1_complex_subset_5_.1_within.csv',
                     header = TRUE, stringsAsFactors = FALSE)

#get 
# list of # of positive GIs per module, 
# list of # of negative GIs per module, 
# purity scores (ranges from -1 to 1, 1 being all positive interaction and -1 being all negative interaction), 
#	purity = actual pos. GI in a module / actual all GI in a module
# list of # of total GIs
observed_pos = data_qGI$No_pos_interactions
observed_neg = data_qGI$No_neg_interactions
observed_purity = data_qGI$Purity_score
prob = sum(observed_pos)/(sum(observed_pos)+sum(observed_neg))

set.seed(10)
#x = rbinom(n = nrow(data_qGI),size = 100,prob = prob)
#background = rescale(x, to = c(-1, 1), from = c(0,100))

# for each module, 
# generate a random positive GI given the # of total GIs and the background probability as weight, from a binomial distribution 
# and calculate the random (background) purity score.
random_pos = c()
random_purity = c()
for(index in c(1:nrow(data_qGI)))
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

# output files

ggsave(file.path('./output_library_3_or_query_1_complex/', 'WB_library_3_or_query_1_complex_subset_5_.1_within_final.pdf'),plot =plot1,height = 14 , width = 14)
