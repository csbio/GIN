# author: Arshia Z. Hassan.

R scripts:
generate_wb_complex_query_vs_library.r
	description:
		generate file with information on within-complex and between-complex GIs.
	input:
		qGI_20211111_nRed.txt
		gi_FDR_20211111_nRed.txt
		data_complex.rda
	output: 
		Within_Between_Enrichment_adjusted_Complex_no_filter_purity_modified.txt
		
create_data_subset_min_interaction.r
	description:
		generate files with subset of within-complex and between-complex GIs based on FDR and No. of interaction tested cutoffs.
	input:
		Within_Between_Enrichment_adjusted_Complex_no_filter_purity_modified.txt
	output: 
		WB_library_3_or_query_1_complex_subset_5_.1_within.csv
		WB_library_3_or_query_1_complex_subset_5_.1_between.csv
		WB_library_3_or_query_1_complex_subset_5_within.csv
		WB_library_3_or_query_1_complex_subset_5_between.csv
		
mitochondrial_gene_fraction_in_complex.R
	description:
		Calculate fraction of mitochondrial genes in each complex.
	input:
		data_complex.rda
		Mitochondial_genelist_1_26_2021_genes.tsv
	output: 
		mitochondrial_fraction_complexes.csv
		
generate_purity_distribution_plot.R
	description:
		Generate distribution plot of purity of negative and positive GIs within complex and between complexes.
	input:
		WB_library_3_or_query_1_complex_subset_5_.1_between.csv
		WB_library_3_or_query_1_complex_subset_5_.1_within.csv
	output:
		WB_library_3_or_query_1_complex_subset_5_.1_between_final.pdf
		WB_library_3_or_query_1_complex_subset_5_.1_within_final.pdf

generate_summary_bar_plot.R
	description:
		Generate barplot of percentage of negative and positive GIs within complex and between complexes.
	input:
		WB_library_3_or_query_1_complex_subset_5_within.csv
		WB_library_3_or_query_1_complex_subset_5_between.csv
	output:
		complex_percentage_bar_plot_wb_5_.pdf
		
generate_purity_distribution_plot_final_wo_mitochondrial_complex_contribution.R
	description:
		Generate distribution of purity of negative and positive GIs within complex and between complexes but removing mitochondrial complexes.
	input:
		WB_library_3_or_query_1_complex_subset_5_.1_between.csv
		WB_library_3_or_query_1_complex_subset_5_.1_within.csv
		mitochondrial_fraction_complexes.csv
	output:
		NO_MITO_.5_WB_library_3_or_query_1_complex_subset_5_.1_between_final.pdf
		NO_MITO_.5_WB_library_3_or_query_1_complex_subset_5_.1_within_final.pdf

generate_barplot_wo_mitochondrial_complex_contribution.R
	description:
		Generate barplot of percentage of negative and positive GIs within complex and between complexes but removing mitochondrial complexes.
	input:
		WB_library_3_or_query_1_complex_subset_5_within.csv
		WB_library_3_or_query_1_complex_subset_5_between.csv
		mitochondrial_fraction_complexes.csv
	output:
		NO_MITO_.5_complex_percentage_bar_plot_wb_5_.pdf
		
generate_complex_overlap_index_version_2_between.R
description:
		calculate overlap index among the complexes to detect redundant complexes.	
	input:
		WB_library_3_or_query_1_complex_subset_5_between.csv
		data_complex.rda
	output:
		subset_5_between_overlap_index_version_2.csv

generate_complex_overlap_index_version_2_within.R
	description:
		calculate overlap index among the complexes to detect redundant complexes.		
	input:
		WB_library_3_or_query_1_complex_subset_2_within.csv
		data_complex.rda
	output:
		subset_2_within_overlap_index_version_2.csv
		
generate_purity_distribution_plot_final_wo_subset_complex_contribution_version_2.R
	description:
		Generate distribution of purity of negative and positive GIs within complex and between complexes but removing redundant complexes.	
	input:
		WB_library_3_or_query_1_complex_subset_5_between.csv
		subset_5_between_overlap_index_version_2.csv
		WB_library_3_or_query_1_complex_subset_5_within.csv
		subset_2_within_overlap_index_version_2.csv
	output:
		NO_SUBSET_.3_WB_library_3_or_query_1_complex_subset_5_.1_between_version_2.pdf
		NO_SUBSET_.3_WB_library_3_or_query_1_complex_subset_5_.1_within_version_2.pdf
		
generate_barplot_wo_subset_complex_contribution_version_2.R
	description:
		Generate barplot of percentage of negative and positive GIs within complex and between complexes but removing redundant complexes.	
	input:
		WB_library_3_or_query_1_complex_subset_5_between.csv
		subset_5_between_overlap_index_version_2.csv
		WB_library_3_or_query_1_complex_subset_5_within.csv
		subset_2_within_overlap_index_version_2.csv
	output:
		NO_SUBSET_.3_complex_percentage_bar_plot_wb_5_version_2.pdf

generate_purity_distribution_plot_final_wo_subset_mito_complex_contribution_version_2.R
	description:
		Generate distribution of purity of negative and positive GIs within complex and between complexes but removing mitochondrial and redundant complexes.	
	input:
		mitochondrial_fraction_complexes.csv
		WB_library_3_or_query_1_complex_subset_5_between.csv
		subset_5_between_overlap_index_version_2.csv
		WB_library_3_or_query_1_complex_subset_5_within.csv
		subset_2_within_overlap_index_version_2.csv
	output:
		NO_MITO_.5_NO_SUBSET_.3_WB_library_3_or_query_1_complex_subset_5_.1_between_version_2.pdf
		NO_MITO_.5_NO_SUBSET_.3_WB_library_3_or_query_1_complex_subset_5_.1_within_version_2.pdf

generate_barplot_wo_subset_mito_complex_contribution_version_2.R
	description:
		Generate barplot of percentage of negative and positive GIs within complex and between complexes but removing mitochondrial and redundant complexes.	
	input:
		mitochondrial_fraction_complexes.csv
		WB_library_3_or_query_1_complex_subset_5_between.csv
		subset_5_between_overlap_index_version_2.csv
		WB_library_3_or_query_1_complex_subset_5_within.csv
		subset_2_within_overlap_index_version_2.csv
	output:
		NO_MITO_.5_NO_SUBSET_.3_complex_percentage_bar_plot_wb_5_version_2.pdf


Run instructions:
	1. Install R and all packages listed at beginning of each script.
	2. Run 'Rscript <script_file_name>' from command line to run each script. Run scripts in the listed order.
	3. Modify file paths in the scripts if necessary.