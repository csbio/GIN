# author: Arshia Z. Hassan. hassa418@umn.edu

local_R_path = "../../local/R/library/"
library("stringi", lib.loc=local_R_path)
library("stringr", lib.loc=local_R_path)
library("rstudioapi", lib.loc=local_R_path)
library("desc", lib.loc=local_R_path)
library("withr", lib.loc=local_R_path)
library("ps", lib.loc=local_R_path)
library("usethis", lib.loc=local_R_path)
library("devtools", lib.loc=local_R_path)
library("RColorBrewer", lib.loc=local_R_path)
library("labeling", lib.loc=local_R_path)
library("farver", lib.loc=local_R_path)
library("digest", lib.loc=local_R_path)
library("ggplot2", lib.loc=local_R_path)
library("ggthemes", lib.loc=local_R_path)
library("ggrepel", lib.loc=local_R_path)

# function to generate PR curve using FLEX 
generate_corum_similarity_plot<-function(corum_sim_list, output_folder, output_file_name_noext, curve_names, cols, title)
{
  
  profile_similarity_plot <- paste(output_file_name_noext, 'PR',  sep = "")
  labs <- c('TP', 'Precision')
  PlotPRSimilarity(corum_sim_list, outfile.name = profile_similarity_plot, 
                   outfile.type = 'pdf', fig.labs = labs, fig.title = title,
                   subsample = T,
                   legend.names = curve_names, 
                   legend.color = cols, save.figure = TRUE)

  files <- c( 
    paste(profile_similarity_plot, '.pdf', sep = "")
  )
  for (f in files)
  {
    file.copy(f, file.path(output_folder, f), overwrite = TRUE)
    file.remove(f)
  }
}

# function to generate AUPRC data using FLEX 
generate_auprc_data<-function(complex_info, data_complex_info, corum_info, output_folder, output_file_name_noext)
{
  complex_info$ID=as.character(complex_info$ID)
  complex_df <- as.data.frame(complex_info, stringsAsFactors = FALSE)
  data_AUPRC <- GetAreaUnderPRCurveForEntities (data_complex_info, corum_info, complex_df)
  auprc_file <- paste('AUPRC_',output_file_name_noext,'.txt', sep = "")
  write.table(data_AUPRC, auprc_file, sep = '\t', row.names = FALSE, quote = FALSE)
  files <- c( 
    auprc_file    
  )
  for (f in files)
  {
    file.copy(f, file.path(output_folder, f), overwrite = TRUE)
    file.remove(f)
  }
}

# load FLEX package
flex_folder <- "./FLEX_R-master"
load_all(flex_folder)
# load CORUM complex standard provided with FLEX
load(file='./FLEX_R-master/data/data_complex.rda')
# generate co-annotation data from CORUM standard
corum <- MakeCoAnnotationFromGeneSymbols(data_standard = data_complex, overlap_length = 1,file_location = "corum")

# load mitochondrial gene list
mito_gene_list_file_genes <- "../../pca_onion/data/Mitochondial_genelist_1_26_2021_genes.tsv"
mito_data_genes <- read.csv(mito_gene_list_file_genes, sep = " ", header = F, stringsAsFactors = F)

# create output directory
output_folder <- './output/'
if (!dir.exists(output_folder)) { dir.create(output_folder, recursive = TRUE) }

###########################
# generate PR-curve plots #
###########################

# create labels and color mappings for PR curve plot
curve_labels <- c('DepMap', 'Gin+DepMap', 'Gin')
col_set <- c()
cols0<-brewer.pal(12, "Paired")
col_set <- c(col_set,cols0[10],cols0[8],cols0[6])

# load depmap data - DepMap version 2020 Q2 pre-processed using the rscript pre_process.R
input_file_name <- "./data/depmap_q2_2020_nona_mean.tsv"
data <- read.csv(input_file_name, sep = "\t", header = TRUE, row.names = 1)

# create similarity network from depmap data
data_t <- data.frame(t(data)) #transpose data - genes as rows
sim_net <- cor(data_t, use = "all.obs",method="pearson") #use pearson correlation to create similarity network

# generate CORUM standard mapping to generate PR curve for depmap data
complex <- CalculatePredictionAndTrueOnLibraryProfiles(corum, sim_net)
corum_sim <- list(list(true = complex$true, predicted = complex$predicted))

# generate AUPRC data against CORUM standard for depmap data and save to file
generate_auprc_data(complex, data_complex, corum, output_folder, 'depmap20q2_complex')

# remove mitochondrial gene pairs from the CORUM standard mapping to generate PR curve for depmap data
complex <- getSubsetOfCoAnnRemovePairs(corum, data.frame(complex), list(mito_data_genes$V1), replace = F)
corum_sim_no_mito <- list(list(true = complex$true, predicted = complex$predicted))

# load bionic integrated features
input_file <- paste("./data/Gin-DepMap-params4_features.tsv", sep='')
data <- read.csv(file.path(input_file), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

# create similarity network bionic integrated features
data_t <- data.frame(t(data))
sim_net <- cor(data_t, use = "all.obs",method="pearson")

# generate CORUM standard mapping to generate PR curve for bionic integrated features
complex <- CalculatePredictionAndTrueOnLibraryProfiles(corum, sim_net) #create mapping for PR curves
corum_sim <- append(corum_sim, list(list(true = complex$true, predicted = complex$predicted)))

# generate AUPRC data against CORUM standard for bionic integrated features and save to file
generate_auprc_data(complex, data_complex, corum, output_folder, 'Gin-DepMap-params4_features_complex')

# remove mitochondrial gene pairs from the CORUM standard mapping to generate PR curve for bionic integrated features
complex <- getSubsetOfCoAnnRemovePairs(corum, data.frame(complex), list(mito_data_genes$V1), replace = F)
corum_sim_no_mito <- append(corum_sim_no_mito, list(list(true = complex$true, predicted = complex$predicted)))

# load GI data
input_file_name <- "./data/qGI_20211111_fullFF.txt"
data <- read.csv(input_file_name, sep = "\t", header = TRUE, row.names = 1,stringsAsFactors=F)

# create similarity network from GI data
data_t <- data.frame(t(data)) #transpose data - genes as rows
sim_net <- cor(data_t, use = "all.obs",method="pearson") #use pearson correlation to create similarity network

# generate CORUM standard mapping to generate PR curve for GI data
complex <- CalculatePredictionAndTrueOnLibraryProfiles(corum, sim_net)
corum_sim <- append(corum_sim, list(list(true = complex$true, predicted = complex$predicted)))

# generate AUPRC data against CORUM standard for GI data and save to file
generate_auprc_data(complex, data_complex, corum, output_folder, 'qGI_20211111_fullFF_complex')

# remove mitochondrial gene pairs from the CORUM standard mapping to generate PR curve for GI data
complex <- getSubsetOfCoAnnRemovePairs(corum, data.frame(complex), list(mito_data_genes$V1), replace = F)
corum_sim_no_mito <- append(corum_sim_no_mito, list(list(true = complex$true, predicted = complex$predicted)))

# generate PR-curve plots

file_ext <- "depmap_gin_bionic-rpco_all gene pairs_CORUM_complex_"
generate_corum_similarity_plot(corum_sim, output_folder, file_ext, curve_labels, col_set, "" )

file_ext <- "depmap_gin_bionic-rpco_no mito. gene pairs_CORUM_complex_"
generate_corum_similarity_plot(corum_sim_no_mito, output_folder, paste(file_ext,sep=''), curve_labels, col_set, "" )

################################################
# prepare data to generate AUPRC scatter plots #
################################################

# load essential gene list
essential_gene_list_file_genes <- "./data/DepMap_essential_20Q2_60_percent.txt"
essential_data_genes <- read.csv(essential_gene_list_file_genes, sep = " ", header = F, stringsAsFactors = F)

orig_folder <- './output/'

# load AUPRC file for depmap
depmap_file <- 'AUPRC_depmap20q2_complex.txt'
AUPRC_depmap <- read.table(file.path(orig_folder,depmap_file), stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
AUPRC_depmap <- AUPRC_depmap[1:4]
colnames(AUPRC_depmap) <- c("ID","Name","Length","AUPRC_depmap")

# load AUPRC file for bionic
gin_file <- 'AUPRC_qGI_20211111_fullFF_complex.txt'
AUPRC_gin <- read.table(file.path(orig_folder,gin_file), stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
AUPRC_gin <- AUPRC_gin[1:4]
colnames(AUPRC_gin) <- c("ID","Name","Length","AUPRC_gin")

# load AUPRC file for GI
int_file <- 'AUPRC_Gin-DepMap-params4_features_complex.txt'
AUPRC_int <- read.table(file.path(orig_folder,int_file), stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
AUPRC_int <- AUPRC_int[1:4]
colnames(AUPRC_int) <- c("ID","Name","Length","AUPRC_bionic")

# merge all AUPRC data
auprc_all <- merge(AUPRC_depmap,AUPRC_gin,all=T)
auprc_all <- merge(auprc_all,AUPRC_int,all=T)

# calculate difference of AUPRC scores between bionic-depmap, bionic-GI, and depmap-GI
auprc_all$subtract_bionic_depmap_auprc = auprc_all$AUPRC_bionic - auprc_all$AUPRC_depmap
auprc_all$subtract_bionic_gin_auprc = auprc_all$AUPRC_bionic - auprc_all$AUPRC_gin
auprc_all$subtract_depmap_gin_auprc = auprc_all$AUPRC_depmap - auprc_all$AUPRC_gin

# for each module in the standard, get the mitochondrial and essential gene overlaps and percentage
id_list <- auprc_all$ID
mito_gene_count = c()
mito_gene_list = c()
essential_gene_list = c()
essential_gene_count = c()
for (cmplx in id_list)
{
  cmplx_gene_list <- data_complex[data_complex$ID == cmplx,]$Genes
  cmplx_gene_list <- unlist(strsplit(cmplx_gene_list, ";"))
  mito_sub <- intersect(cmplx_gene_list,mito_data_genes$V1)
  mito_sub_str <- paste(mito_sub,collapse=";")
  mito_gene_list = c(mito_gene_list,mito_sub_str)
  mito_gene_count = c(mito_gene_count,length(mito_sub)/length(cmplx_gene_list))
  essential_sub <- intersect(cmplx_gene_list,essential_data_genes$V1)
  essential_sub_str <- paste(essential_sub,collapse=";")
  essential_gene_list = c(essential_gene_list,essential_sub_str)
  essential_gene_count = c(essential_gene_count,length(essential_sub)/length(cmplx_gene_list))
}
auprc_all$Fraction_mito =   mito_gene_count
auprc_all$Mito_genes = mito_gene_list
auprc_all$Fraction_essential =   essential_gene_count
auprc_all$Essential_genes = essential_gene_list

# save merged AUPRC data
write.csv(auprc_all,"./output/AUPRC_merged_depmap_gin_bionic_rpco_complex.csv",quote = T)

#############################################
# generate scatter plot -- bionic vs depmap #
############################################# 

#auprc_all = read.csv("./output/AUPRC_merged_depmap_gin_bionic_rpco_complex.csv")
#output_folder <- './output/'
data = auprc_all
plot_df = data[, c('AUPRC_depmap', 'AUPRC_bionic')] 
plot_df <- plot_df[complete.cases(plot_df), ]
diff_cut_off = .2
auprc_cut_off = .35
plot_df_1 = data[which((abs(data$subtract_bionic_depmap_auprc)> diff_cut_off)
                       & ((data$AUPRC_depmap > auprc_cut_off) | (data$AUPRC_bionic > auprc_cut_off))
), c('AUPRC_depmap', 'AUPRC_bionic', 'Name', 'ID')]
plot_df_1 <- plot_df_1[complete.cases(plot_df_1), ]
# selected = c('CDK8 subcomplex (CCNC, CDK8, MED12, MED13)',
#              'DNA ligase IV-XRCC4-XLF complex',
#              'KICSTOR complex',
#              'CBF-DNA complex',
#              'HUSH complex',
#              'PBAF complex (Polybromo- and BAF containing complex)',
#              'Tetrameric COG subcomplex',
#              'BORC complex',
#              'Cofilin-actin-CAP1 complex',
#              'EARP complex',
#              'TSC complex'
# )#complex
selected = c('7583',
             '4478',
             '359',
             '6754',
             '6891',
             '556',
             '3525',
             '7545',
             '2255',
             '6459',
             '6889'
)#complex
plot_df_2 <- plot_df_1[which(plot_df_1$ID %in% selected),]

plot_main <- ggplot(data = plot_df,
                    aes(x=AUPRC_depmap, y=AUPRC_bionic)) +
  geom_point(size=4, alpha = 1, color = 'grey') +
  geom_point(data = plot_df_1, size=4, alpha = 1, 
             aes(x=AUPRC_depmap, y=AUPRC_bionic)) +
  geom_abline(intercept = diff_cut_off, slope = 1, linetype = "dashed") +
  geom_abline(intercept = -diff_cut_off, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = auprc_cut_off, linetype = "dashed") +
  geom_vline(xintercept = auprc_cut_off, linetype = "dashed") +
  geom_point(data = plot_df_2, size=4, alpha = 1, color = 'red',
             aes(x=AUPRC_depmap, y=AUPRC_bionic)) +
  geom_label_repel(data = plot_df_2, max.overlaps = 50, alpha=.3,
                   aes(
                     x=AUPRC_depmap, y=AUPRC_bionic,label = Name, size = NULL, color = NULL), 
                   size = 5

  )+
  xlab("DepMap") +
  ylab("Bionic(GI,DepMap)") +
  ylim(0, 1) +
  xlim(0, 1) +
  theme_bw() +
  labs(colour = "Domain") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size = 2),
        legend.title = element_text(size = 20,colour = "black"),
        legend.text = element_text(size = 20,colour = "black"),
        axis.title = element_text(size = 30,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text = element_text(size = 30,colour = "black")
        ,legend.position = "none"
  )

ggsave(file.path(output_folder, "scatter_auprc_depmap_bionic_complex.pdf"),height = 10 , width = 10 )

##########################################
# generate scatter plot --  bionic vs gi # 
##########################################

data = auprc_all
plot_df = data[, c('AUPRC_gin', 'AUPRC_bionic')]  
plot_df <- plot_df[complete.cases(plot_df), ]
diff_cut_off = .2
auprc_cut_off = .35
plot_df_1 = data[which((abs(data$subtract_bionic_gin_auprc)> diff_cut_off)
                        & ((data$AUPRC_gin > auprc_cut_off) | (data$AUPRC_bionic > auprc_cut_off))
), c('AUPRC_gin', 'AUPRC_bionic', 'Name')]
plot_df_1 <- plot_df_1[complete.cases(plot_df_1), ]

plot_main <- ggplot(data = plot_df,
                    aes(x=AUPRC_gin, y=AUPRC_bionic)) +
  geom_point(size=4, alpha = 1, color = 'grey') +
  geom_point(data = plot_df_1, size=4, alpha = 1,
             aes(x=AUPRC_gin, y=AUPRC_bionic)) +
  geom_abline(intercept = diff_cut_off, slope = 1, linetype = "dashed") +
  geom_abline(intercept = -diff_cut_off, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = auprc_cut_off, linetype = "dashed") +
  geom_vline(xintercept = auprc_cut_off, linetype = "dashed") +
  xlab("GI") +
  ylab("Bionic(GI,DepMap)") +
  ylim(0, 1) +
  xlim(0, 1) +
  theme_bw() +
  labs(colour = "Domain") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size = 2),
        legend.title = element_text(size = 20,colour = "black"),
        legend.text = element_text(size = 20,colour = "black"),
        axis.title = element_text(size = 30,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text = element_text(size = 30,colour = "black")
        ,legend.position = "none"
  )

ggsave(file.path(output_folder, "scatter_auprc_gin_bionic_complex.pdf"),height = 10 , width = 10 )

##########################################
# generate scatter plot --  gi vs depmap #
##########################################

data = auprc_all
plot_df = data[, c('AUPRC_depmap', 'AUPRC_gin')]
plot_df <- plot_df[complete.cases(plot_df), ]
diff_cut_off = .2
auprc_cut_off = .35
plot_df_1 = data[which((abs(data$subtract_depmap_gin_auprc)> diff_cut_off)
                       & ((data$AUPRC_depmap > auprc_cut_off) | (data$AUPRC_gin > auprc_cut_off))
), c('AUPRC_depmap', 'AUPRC_gin', 'Name')]
plot_df_1 <- plot_df_1[complete.cases(plot_df_1), ]

plot_main <- ggplot(data = plot_df,
                    aes(x=AUPRC_depmap, y=AUPRC_gin)) +
  geom_point(size=4, alpha = 1, color = 'grey') +
  geom_point(data = plot_df_1, size=4, alpha = 1,
             aes(x=AUPRC_depmap, y=AUPRC_gin)) +
  geom_abline(intercept = diff_cut_off, slope = 1, linetype = "dashed") +
  geom_abline(intercept = -diff_cut_off, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = auprc_cut_off, linetype = "dashed") +
  geom_vline(xintercept = auprc_cut_off, linetype = "dashed") +
  ylab("GI") +
  xlab("DepMap") +
  ylim(0, 1) +
  xlim(0, 1) +
  theme_bw() +
  labs(colour = "Domain") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size = 2),
        legend.title = element_text(size = 20,colour = "black"),
        legend.text = element_text(size = 20,colour = "black"),
        axis.title = element_text(size = 30,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text = element_text(size = 30,colour = "black")
        ,legend.position = "none"
  )

ggsave(file.path(output_folder, "scatter_auprc_depmap_gin_complex.pdf"),height = 10 , width = 10 )
