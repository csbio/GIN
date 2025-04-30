# author: Arshia Z. Hassan.

R scripts:

generate_siginificant_GIs_files.R
	description: Generate files with significant positive, negative GIs, and all GI gene pair list.
	input: 
		qGI_20211111_nRed.txt
		gi_FDR_20211111_nRed.txt
	output:
		all_significant_positive_GIs_gene_pairs.csv
		all_significant_negative_GIs_gene_pairs.csv
		all_GIs_gene_pairs.csv
		
percent_GI_vs_sequence_id_ensemble_set.R
	description: generate line plots - percent of positive and negative GI vs sequence ID bin for various expression filters from ensemble paralog set.
	input: 
		ensembl_ohnolog_pairs_complete.rds
		Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx
		all_significant_positive_GIs_gene_pairs.csv
		all_significant_negative_GIs_gene_pairs.csv
		all_GIs_gene_pairs.csv
	output:
		no_exp_vs_both_exp_vs_one_exp__exp_id_open_bin_dot_plot_.pdf
		
percent_GI_vs_recalculated_family_size_ensemble_set.R
	description: Generate plot: percentage of negative GI vs family size from ensemble paralog set (sequence ID > 50).
	input: 
		ensembl_ohnolog_pairs_complete.rds
		Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx
		all_significant_negative_GIs_gene_pairs.csv
		all_GIs_gene_pairs.csv
	output:
		barplot_id_open_bin_50_neg_0.pdf
		
percent_GI_vs_recalculated_family_size_ohnolog_set.R
	description: Generate plot: percentage of negative GI vs family size from ohnolog paralog set.
	input: 
		hsapiens.Pairs.Relaxed.2R.txt
		Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx
		all_significant_negative_GIs_gene_pairs.csv
		all_GIs_gene_pairs.csv
	output:
		barplot_both_exp_open_bin_neg.pdf

Run instructions:
	1. Install R and all packages listed at beginning of each script.
	2. Run 'Rscript <script_file_name>' from command line to run each script.