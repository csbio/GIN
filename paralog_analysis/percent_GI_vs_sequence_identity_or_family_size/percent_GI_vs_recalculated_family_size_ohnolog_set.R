library("withr")
library("farver")
library("labeling")
library("digest")
library("readxl")
library("dplyr")
library("reshape2")
library("ggplot2")
library("ggthemes")
library("gridExtra")
library("forcats")
library("igraph")

output_folder = "./"

# load significant negative GIs and all screened GIs
all_significant_negative_GIs = read.csv("./all_significant_negative_GIs_gene_pairs.csv")
all_GIs = read.csv("./all_GIs_gene_pairs.csv")

# calculate percentage of significant negative GIs
background_neg = nrow(all_significant_negative_GIs)/nrow(all_GIs)*100 

# load ohnolog gene pair list
ohnolog = read.csv('./hsapiens.Pairs.Relaxed.2R.txt',sep='\t')

ohnolog = ohnolog %>%
  rowwise() %>%
  mutate(sorted_pair = paste(unlist(sort(c(as.character(Symbol1),Symbol2))), collapse = "_"))

# filter paralog gene pair list based on expression - keep gene pairs with only log2(TPM+1) expression > 1
supp_table_2_file = "./Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx"
supp_table_2 = read_excel(supp_table_2_file)
supp_table_sub = supp_table_2[which(supp_table_2$`HAP1_mRNA expression [log2(TPM+1)]`>1),]
ohnolog_sub_ = ohnolog[which(ohnolog$Symbol1 %in% supp_table_sub$`Gene Symbol` & ohnolog$Symbol2 %in% supp_table_sub$`Gene Symbol`),]

# recalculate family sizes with respect to the filtered paralog set 
sub_ad = as.matrix(get.adjacency(graph.data.frame(ohnolog_sub_[,c('Symbol1','Symbol2')])))
graf  <- as.undirected(graph.adjacency(sub_ad))
clusters = groups(components(graf))

fam_size_df = data.frame(matrix(nrow=0,ncol=2))

for(i in 1:length(clusters))
{
  clus = clusters[[i]]
  fam_size = length(clus)
  for(g in clus)
  {
    fam_size_df = rbind(fam_size_df,c(g,fam_size))
  }
}
colnames(fam_size_df) = c('gene','family_size')

ohnolog_sub_$family_size = 0
for(i in 1:nrow(ohnolog_sub_))
{
  ohnolog_sub_[i,'family_size'] = as.numeric(fam_size_df$family_size[which(fam_size_df$gene %in% ohnolog_sub_[i,'Symbol1'])])
}

# filter based on screened GI pairs
ohnolog_sub_ = ohnolog_sub_[which(ohnolog_sub_$sorted_pair %in% all_GIs$x),]
ohnolog_sub_ = ohnolog_sub_[!duplicated(ohnolog_sub_$sorted_pair),]

# get distribution of family sizes
ohnolog_sub_freq = data.frame(table(ohnolog_sub_$family_size))

ohnolog_sub_neg = ohnolog_sub_[which(ohnolog_sub_$sorted_pair %in% all_significant_negative_GIs$x),]


###generate plot data collapsed################################################################################
family_size_list = as.numeric(intersect(c(2,3,4),ohnolog_sub_freq$Var1))

fraction_list = c()
total_list =c()

# get negative GI percentage with family size of 2,3,4
for(family in family_size_list)
{
  family = as.numeric(family)
  neg = sum(ohnolog_sub_neg$family_size==family)/sum(ohnolog_sub_$family_size==family)*100
  fraction_list = c(fraction_list,neg)
  total_list = c(total_list, sum(ohnolog_sub_$family_size==family))
}

# get negative GI percentage with family size of greater than 4
if(any(as.numeric(levels(ohnolog_sub_freq$Var1))[ohnolog_sub_freq$Var1] >4)){
  family_size_list = c(family_size_list, '>4')
  neg = sum(ohnolog_sub_neg$family_size >4)/sum(ohnolog_sub_$family_size >4)*100
  fraction_list = c(fraction_list,neg)
  total_list = c(total_list, sum(ohnolog_sub_$family_size >4))
}

# get negative GI percentage with family size of greater than 0 
family_size_list = c(family_size_list, 0)
neg = nrow(ohnolog_sub_neg)/nrow(ohnolog_sub_)*100
fraction_list = c(fraction_list,neg)
total_list = c(total_list, nrow(ohnolog_sub_))

df_plot = data.frame(
  'family_size' = family_size_list,
  'total' = total_list,
  'fraction' = fraction_list,
  'GI' = rep(c( "neg" ) , length(family_size_list))
)


df_plot$x_label = paste(df_plot$family_size,' (',df_plot$total,')',sep='')

# generate plot percentage of GI vs family size

PVEplot <- ggplot(data = df_plot, aes(fill=GI, y=fraction, x=fct_inorder(x_label))) +
  geom_bar(position="dodge", stat="identity") +
  geom_hline(yintercept = background_neg, color='blue',linetype=2) +
  xlab("Family size") +
  ylab("% of GI") +
  theme_bw() +
  scale_fill_manual(values=c('blue','yellow')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), #axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        legend.title = element_text(family="Sans",size = 50,colour = "black"),
        legend.text = element_text(family="Sans",size = 50,colour = "black"),
        axis.title = element_text(family="Sans",size = 30,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(family="Sans",size = 30,colour = "black"),
        axis.text.x = element_text(family="Sans",size = 30,colour = "black",angle = 90, vjust = 0.5),
        plot.margin = unit(c(1,1,1,1), "cm"),
        axis.ticks.length=unit(0.1,"inch")
        #, legend.position="none"
  )

ggsave(file.path(output_folder, paste0("barplot_both_exp_open_bin_neg.pdf")), device = cairo_pdf,height = 10 , width = 20 )

