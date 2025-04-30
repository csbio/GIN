# load CORUM complex information file
load("data_complex.rda")

####within complex################################################################################

#load within complex analysis output file: WB_library_3_or_query_1_complex_subset_5_within_.csv 
data_qGI <- read.csv('output_library_3_or_query_1_complex/WB_library_3_or_query_1_complex_subset_2_within.csv',
                     header = TRUE, stringsAsFactors = FALSE)

# for each complex in the within complex analysis output file, get the complex size in terms of # of member genes
# required to get a list of complexes ordered by size
comlex_list <- unique(data_qGI$Column_Complex)
complex_df = data.frame(matrix(nrow=0,ncol=2))
colnames(complex_df) = c('complex','length')
for(i in 1 : (length(comlex_list)) )
{
  gene_list_i_ = unlist(data_complex$Genes[which(data_complex$Name==comlex_list[i])]) #get gene list from CORUM complex information file
  gene_list_i = unlist(strsplit(gene_list_i_, ';')) # get gene list from ; separated string
  gene_list_i = gsub(' ', '', gene_list_i) # Replacing any spaces with nothing
  gene_list_i = sort(unique(gene_list_i)) #Remove duplicate genes
  temp_df = data.frame(complex=c(comlex_list[i]), length=c(length(gene_list_i)))
  complex_df = rbind(complex_df,temp_df)
}
# order complexes by size
complex_df = complex_df[order(complex_df$length,decreasing = T),]
# get complex list ordered by size
comlex_list = complex_df$complex

# create dataframe to store overlap index per complex pair
sim_jacrd = data.frame(matrix(nrow=0,ncol=8))
colnames(sim_jacrd) = c('com_1','com_2','overlap_index','length_1','length_2','length_diff','gene_1','gene_2')
# for each complex from list ordered by size, calculate overlap index per complex pair
for (i in 1 : (length(comlex_list) - 1) )
{
  gene_list_i_ = unlist(data_complex$Genes[which(data_complex$Name==comlex_list[i])]) #get gene list from CORUM complex information file
  gene_list_i = unlist(strsplit(gene_list_i_, ';')) # get gene list from ; separated string
  gene_list_i = gsub(' ', '', gene_list_i) # Replacing any spaces with nothing
  gene_list_i = sort(unique(gene_list_i)) #Remove duplicate genes
  gene_list_i_ = paste0(gene_list_i, collapse = ';') # generate gene string separated by ;
  for (j in (i+1):length(comlex_list))
  {
    gene_list_j_= data_complex$Genes[which(data_complex$Name==comlex_list[j])] #get gene list from CORUM complex information file
    gene_list_j = unlist(strsplit(gene_list_j_, ';')) # get gene list from ; separated string
    gene_list_j = gsub(' ', '', gene_list_j) # Replacing any spaces with nothing
	gene_list_j = sort(unique(gene_list_j)) #Remove duplicate genes
    gene_list_j_ = paste0(gene_list_j, collapse = ';') # generate gene string separated by ;
    # Overlap Index: (size of intersection of two complexes) / (minimum of the two complex sizes)
    jaccard <- length(intersect(gene_list_i, gene_list_j)) / min(length(gene_list_i), length(gene_list_j))
	# Append overlap index to dataframe
    temp_sim_jacrd = data.frame(com_1=comlex_list[i], com_2=comlex_list[j], overlap_index = c(jaccard), 
                                length_1=c(length(gene_list_i)),length_2=c(length(gene_list_j)),
                                length_diff=c(abs(length(gene_list_i)-length(gene_list_j))),
                                gene_1 = c(gene_list_i_), gene_2=c(gene_list_j_))
    sim_jacrd = rbind(sim_jacrd,temp_sim_jacrd)
  }
}

#save to file
write.csv(sim_jacrd,'output_library_3_or_query_1_complex/subset_2_within_overlap_index_version_2.csv',row.names = F)



