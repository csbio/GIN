
#load mitochondrial gene list
mito_genes <- read.table('./Mitochondial_genelist_1_26_2021_genes.tsv')

# load CORUM complex information file
load("./data_complex.rda")

comlex_list <- data_complex$Name

# create dataframe to store mitochondrial gene fraction, count, and gene names per complex 
sim_jacrd = data.frame(matrix(nrow=0,ncol=6))
colnames(sim_jacrd) = c('complex','fraction_mito','length','length_mito','gene','gene_mito')
for (i in 1 : (length(comlex_list)) )
{
  gene_list_i_ = unlist(data_complex$Genes[which(data_complex$Name==comlex_list[i])]) #get gene list from CORUM complex information file
  gene_list_i = unlist(strsplit(gene_list_i_, ';')) # get gene list from ; separated string
  gene_list_i = gsub(' ', '', gene_list_i) # Replacing any spaces with nothing
  gene_list_i_ = paste0(gene_list_i, collapse = ';') # generate gene string separated by ;
  
  gene_list_mito = intersect(gene_list_i, mito_genes$V1) # get mitochondrial genes in the complex
  gene_list_mito_ = paste0(gene_list_mito, collapse = ';') # generate gene string separated by ;
  # Fraction of mitochondrial genes in the complex
  fraction <-  length(gene_list_mito) / length(gene_list_i)

  temp_sim_jacrd = data.frame(complex=c(comlex_list[i]), fraction_mito = c(fraction), 
                              length=c(length(gene_list_i)),length_mito=c(length(gene_list_mito)),
                              gene = c(gene_list_i_), gene_mito=c(gene_list_mito_))
  sim_jacrd = rbind(sim_jacrd,temp_sim_jacrd)

}

# discard any complexes with no mitochondrial genes before saving to file
sim_jacrd = sim_jacrd[which(!sim_jacrd$length_mito==0),]
# save output to file
write.csv(sim_jacrd,'./output_library_3_or_query_1_complex/mitochondrial_fraction_complexes.csv',row.names = F)
