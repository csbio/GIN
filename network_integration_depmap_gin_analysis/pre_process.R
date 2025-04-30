# author: Arshia Z. Hassan. hassa418@umn.edu

start_time <- Sys.time()
split_gene_name<-function(x)
{
  x1<-strsplit(x,"\\.")[[1]][1]
  return( x1)
}

# load depmap data: version 2020 Q2
input_file_name <- "./data/Achilles_gene_effect.csv"
Achilles_gene_effect <- read.csv(file = input_file_name,row.names = 1, header = TRUE)

# extract gene names from column names and update column names with gene names
cols<-c(colnames(Achilles_gene_effect))
cols2<-unlist(lapply(cols,split_gene_name))
colnames(Achilles_gene_effect) <- c(cols2)

# transpose data cell-lines X gene to gene X cell-lines
Achilles_gene_effect_t <- data.frame(t(Achilles_gene_effect))

# replace NA values with per-gene average
k <- which(is.na(Achilles_gene_effect_t), arr.ind=TRUE)
Achilles_gene_effect_t[k] <- rowMeans(Achilles_gene_effect_t, na.rm=TRUE)[k[,1]]# replace with row average

# save data to file 
write.table(Achilles_gene_effect_t,file="./data/depmap_q2_2020_nona_mean.tsv",sep='\t',row.names = TRUE,col.names = TRUE, quote=FALSE)

end_time <- Sys.time()
print("Time: ")
print(end_time - start_time)

print("Memory profile: ")
print(memory.profile())

print("Object-wise Memory profile: ")

sapply(ls(), function(x) {print(object.size(get(x)),standard = "legacy", units="Gb") })

