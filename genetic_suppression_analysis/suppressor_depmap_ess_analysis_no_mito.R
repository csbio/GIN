# author: Arshia Z. Hassan. hassa418@umn.edu
options(java.parameters = "- Xmx1024m")
#local_R_path = "../../../local/R/library/"
library("withr")
library("readxl")
library("ggthemes")
library("ggplot2")
library("reshape2")
library("grid")
library("ggpubr")
library("rstatix")

gene_info = read_excel("./Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx")

#load all significant positive GIs and suppression score
supp_file = "./all_significant_positive_GIs_supp_metric_2_.9_.csv"
supp_data = read.csv(supp_file, stringsAsFactors=FALSE)

mito_list = read.csv('./Mitochondial_genelist_1_26_2021_genes.tsv',sep = '\t',stringsAsFactors = F,header=F)

supp_cutoff_list = c(.2,.5,.8)
for (supp_cutoff in supp_cutoff_list)
{
	# define suppressors based on cutoff
	supp_degree = data.frame(table(supp_data$library_gene[which(supp_data$`score`>supp_cutoff)]))
	all_pos_degree = data.frame(table(supp_data$library_gene))

	#all Hap1 essential genes
	all_hap1_ess_gene = gene_info$`Gene Symbol`[which(gene_info$`HAP1_essential gene`==1 )]
	all_hap1_ess_gene = setdiff(all_hap1_ess_gene,mito_list$V1)
	all_hap1_ess_gene_score = gene_info[which(gene_info$`Gene Symbol` %in% all_hap1_ess_gene),c('Gene Symbol',"fractionEssDepMap")]
	mean(all_hap1_ess_gene_score$fractionEssDepMap,na.rm=T)
	median(all_hap1_ess_gene_score$fractionEssDepMap,na.rm=T)

	#all Hap1 essential positive GI with at least one interaction 
	hap1_ess_pos_deg_1_gene = intersect(all_hap1_ess_gene,all_pos_degree$Var1[which(all_pos_degree$Freq>0)])
	hap1_ess_pos_deg_1_gene = setdiff(hap1_ess_pos_deg_1_gene,mito_list$V1)
	hap1_ess_pos_deg_1_gene_score = gene_info[which(gene_info$`Gene Symbol` %in% hap1_ess_pos_deg_1_gene),c('Gene Symbol',"fractionEssDepMap")]
	mean(hap1_ess_pos_deg_1_gene_score$fractionEssDepMap,na.rm=T)
	median(hap1_ess_pos_deg_1_gene_score$fractionEssDepMap,na.rm=T)

	#all Hap1 essential suppressor with at least one interaction 
	hap1_ess_supp_deg_0_gene = intersect(all_hap1_ess_gene,supp_degree$Var1[which(supp_degree$Freq >0)])
	hap1_ess_supp_deg_0_gene = setdiff(hap1_ess_supp_deg_0_gene,mito_list$V1)
	hap1_ess_supp_deg_0_gene_score = gene_info[which(gene_info$`Gene Symbol` %in% hap1_ess_supp_deg_0_gene),c('Gene Symbol',"fractionEssDepMap")]
	mean(hap1_ess_supp_deg_0_gene_score$fractionEssDepMap,na.rm=T)
	median(hap1_ess_supp_deg_0_gene_score$fractionEssDepMap,na.rm=T)

	# Hap1 essential non-suppressor
	hap1_ess_pos_deg_0_gene = setdiff(all_hap1_ess_gene,mito_list$V1)
	hap1_ess_non_supp_gene = setdiff(hap1_ess_pos_deg_0_gene,hap1_ess_supp_deg_0_gene)
	hap1_ess_non_supp_gene = setdiff(hap1_ess_non_supp_gene,mito_list$V1)
	hap1_ess_non_supp_gene_score = gene_info[which(gene_info$`Gene Symbol` %in% hap1_ess_non_supp_gene),c('Gene Symbol',"fractionEssDepMap")]
	mean(hap1_ess_non_supp_gene_score$fractionEssDepMap,na.rm=T)
	median(hap1_ess_non_supp_gene_score$fractionEssDepMap,na.rm=T)

	wt = merge(all_hap1_ess_gene_score,hap1_ess_pos_deg_1_gene_score, by='Gene Symbol', all=T)
	colnames(wt) = c('Gene Symbol', 'Hap1_ess.','Hap1_ess._pos_GI')
	wt = merge(wt,hap1_ess_supp_deg_0_gene_score, by='Gene Symbol', all=T)
	colnames(wt) = c('Gene Symbol', 'Hap1_ess.','Hap1_ess._pos_GI', 'Hap1_ess._supp.')
	wt = merge(wt,hap1_ess_non_supp_gene_score, by='Gene Symbol', all=T)
	colnames(wt) = c('Gene Symbol', 'Hap1_ess.','Hap1_ess._pos_GI', 'Hap1_ess._supp.', 'Hap1_ess._non_supp.')

	wt_all = melt(wt[,c('Hap1_ess.','Hap1_ess._pos_GI', 'Hap1_ess._supp.', 'Hap1_ess._non_supp.')])
	wt_all =wt_all[complete.cases(wt_all),]
	colnames(wt_all) = c('subset','score')

	stat.test <- wt_all %>%
	  wilcox_test(score~subset, p.adjust.method = "none",detailed=T) 
	stat.test <- stat.test %>%
	  add_xy_position()

	count_df = data.frame(table(wt_all$subset))

	PVEplot <- ggplot(data = wt_all, aes(x = factor(subset), y = score)) +
	  geom_boxplot() +
	  xlab("") +
	  ylab("Fraction of essential DepMap") +
	  scale_x_discrete(labels= c('Hap1 ess.','Hap1 ess. pos. GI degree>0', 'Hap1_ess. supp.', 'Hap1 ess. non-supp.')) +
	  stat_pvalue_manual(
		stat.test,  label = "{p}", tip.length = 0.02, hide.ns = TRUE, #"{p}{p.adj.signif}"
		step.increase = 0.08, bracket.nudge.y = 0.08, size=10
	  )+
	  geom_text(data = count_df, aes(x=c(1,2,3,4), y=-.1,label = Freq),  size=10, position = position_dodge(0.9), vjust = 0)+
	  theme_bw() +
	  scale_colour_brewer(palette = "Set1") +
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
			, legend.position="none"
	  )+ 
	  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

	  
	output_folder = "./"
	ggsave(file.path(output_folder, paste0("boxplot_fraction_of_esssential_depmap_suppressor_no_mito_", supp_cutoff, ".pdf")), device = cairo_pdf,height = 15 , width = 10 )
}