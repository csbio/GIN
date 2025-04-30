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


#filter based on expression
supp_table_2_file = "./Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx"
supp_table_2 = read_excel(supp_table_2_file)
supp_table_sub = supp_table_2[which(supp_table_2$`HAP1_mRNA expression [log2(TPM+1)]`>1),]

#######################################################################
degree_data = read.csv('./qGI_degree_auc_ess.txt',sep = '\t')
density_data = degree_data[,c("gene","qGI_degree_std","qGI_degree_std_neg","qGI_degree_std_pos")]
density_data$qGI_degree_std = density_data$qGI_degree_std/222 #divide by the total # of non-redundant query genes
density_data$qGI_degree_std_neg = density_data$qGI_degree_std_neg/222
density_data$qGI_degree_std_pos = density_data$qGI_degree_std_pos/222

#######################################################################

ohnolog = readRDS(file='./ensembl_ohnolog_pairs_complete.rds')
ohnolog = ohnolog[!duplicated(ohnolog$sorted_pair),]
ohnolog_sub_ = ohnolog[which(ohnolog$name1 %in% supp_table_sub$`Gene Symbol` & ohnolog$name2 %in% supp_table_sub$`Gene Symbol`),]

#filter based on ident_min
iden_cutoff = 20#50
ohnolog_sub_ = ohnolog_sub_[which(ohnolog_sub_$ident_max > iden_cutoff),]

#recalculate family sizes
sub_ad = as.matrix(get.adjacency(graph.data.frame(ohnolog_sub_[,c('name1','name2')])))
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

for(i in 1:nrow(ohnolog_sub_))
{
  ohnolog_sub_[i,'family_size'] = as.numeric(fam_size_df$family_size[which(fam_size_df$gene == ohnolog_sub_[i,'name1'])])
}

ohnolog_sub_ = ohnolog_sub_[,c('name1','name2','family_size')]
ohnolog_sub_2 = ohnolog_sub_[,c('name1','family_size')]
ohnolog_sub_3 = ohnolog_sub_[,c('name2','family_size')]
colnames(ohnolog_sub_3) = c('name1','family_size')
ohnolog_sub_2 = rbind(ohnolog_sub_2,ohnolog_sub_3)
ohnolog_sub_2 = unique(ohnolog_sub_2)


#filter based on screened genes
ohnolog_sub_2 = ohnolog_sub_2[which(ohnolog_sub_2$name1 %in% density_data$gene),]

ohnolog_sub_2$qGI_degree_std_pos = density_data$qGI_degree_std_pos[match(ohnolog_sub_2$name1,density_data$gene)]
ohnolog_sub_2$qGI_degree_std_neg = density_data$qGI_degree_std_neg[match(ohnolog_sub_2$name1,density_data$gene)]

ohnolog_sub_data = ohnolog_sub_2 %>% group_by(family_size) %>% summarise(qGI_degree_std_pos= mean(qGI_degree_std_pos),qGI_degree_std_neg= mean(qGI_degree_std_neg))


###generate plot data collapsed################################################################################
family_size_list = as.numeric(intersect(c(2,3,4),ohnolog_sub_2$family_size))

fraction_list = c()
total_list = c()
for(family in family_size_list)
{
  family = as.numeric(family)
  pos = mean(ohnolog_sub_2$qGI_degree_std_pos[ohnolog_sub_2$family_size==family])*100
  neg = mean(ohnolog_sub_2$qGI_degree_std_neg[ohnolog_sub_2$family_size==family])*100
  fraction_list = c(fraction_list,pos,neg)
  total_list = c(total_list, sum(ohnolog_sub_2$family_size==family))
}

if(any(ohnolog_sub_2$family_size >4)){
  family_size_list = c(family_size_list, '>4')
  pos = mean(ohnolog_sub_2$qGI_degree_std_pos[ohnolog_sub_2$family_size >4])*100
  neg = mean(ohnolog_sub_2$qGI_degree_std_neg[ohnolog_sub_2$family_size >4])*100
  fraction_list = c(fraction_list,pos,neg)
  total_list = c(total_list, sum(ohnolog_sub_2$family_size >4))
}

family_size_list = c(family_size_list, 0)
pos = mean(ohnolog_sub_2$qGI_degree_std_pos)*100
neg = mean(ohnolog_sub_2$qGI_degree_std_neg)*100
fraction_list = c(fraction_list,pos,neg)
total_list = c(total_list, nrow(ohnolog_sub_2))

df_plot = data.frame(
  'family_size' = rep(family_size_list,each = 2),
  'fraction' = fraction_list,
  'GI' = rep(c("pos" , "neg" ) , length(family_size_list)),
  'total' = rep(total_list,each = 2)
)


df_plot$x_label = paste(df_plot$family_size,' (',df_plot$total,')',sep='')


###generate plot#################################################################################


PVEplot <- ggplot(data = df_plot, aes(fill=GI, y=fraction, x=fct_inorder(x_label))) +
  geom_bar(position="dodge", stat="identity") +
  xlab("Family size") +
  ylab("Average density") +

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

ggsave(file.path(output_folder, paste0("ensemble_degree_density_vs_family_size_barplot_id_open_bin_", iden_cutoff, ".pdf")), device = cairo_pdf,height = 10 , width = 20 )
