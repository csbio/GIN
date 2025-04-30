
data_WB <- read.table('./output_library_3_or_query_1_complex/Within_Between_Enrichment_adjusted_Complex_no_filter_purity_modified.txt',
                       header = TRUE, sep = '\t', stringsAsFactors = FALSE)
					   
data_WB_within <- data_WB[which(data_WB$type=='within' & data_WB$Enrichment_total<.1 & data_WB$No_of_interaction_tested>5),]
write.csv(data_WB_within,file.path('./output_library_3_or_query_1_complex/', "WB_library_3_or_query_1_complex_subset_5_.1_within.csv"),row.names = F)

data_WB_between <- data_WB[which(data_WB$type=='between' & data_WB$Enrichment_total<.1 & data_WB$No_of_interaction_tested>5),]
write.csv(data_WB_between,file.path('./output_library_3_or_query_1_complex/', "WB_library_3_or_query_1_complex_subset_5_.1_between.csv"),row.names = F)

data_WB_within <- data_WB[which(data_WB$type=='within' & data_WB$No_of_interaction_tested>5),]
write.csv(data_WB_within,file.path('./output_library_3_or_query_1_complex/', "WB_library_3_or_query_1_complex_subset_5_within.csv"),row.names = F)

data_WB_between <- data_WB[which(data_WB$type=='between' & data_WB$No_of_interaction_tested>5),]
write.csv(data_WB_between,file.path('./output_library_3_or_query_1_complex/', "WB_library_3_or_query_1_complex_subset_5_between.csv"),row.names = F)

