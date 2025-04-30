# author: Arshia Z. Hassan.

Asymmetry analysis R scripts:

generate_asymmetry_distribution_plot_ohnolog_set.R
	description: Generate distribution plot: ratio of negative interaction degree of paralog gene pairs from ohnolog set
	input: 
		hsapiens.Pairs.Relaxed.2R.txt
		Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx
		qGI_degree_auc_ess.txt
	output:
		density_negative_GI_ratio_ohnolog_hsapiens.Pairs.Relaxed.2R_.pdf
		
generate_asymmetry_distribution_plot_ensemble_set.R
	description: Generate distribution plot: ratio of negative interaction degree of paralog gene pairs (sequence id >50) from ensemble set
	input: 
		ensembl_ohnolog_pairs_complete.rds
		Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx
		qGI_degree_auc_ess.txt
	output:
		density_negative_GI_ratio_seqid_50_.pdf
		
generate_asymmetry_null_model_empirical_pval_ensemble_set.R
	description: Generate null model and calculate empirical p-val comparing null model and observed ratios: ratio of negative interaction degree of paralog gene pairs (sequence id >50) from ensemble set
	input: 
		hsapiens.Pairs.Relaxed.2R.txt
		Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx
		qGI_degree_auc_ess.txt
	output:
		average_null_negative_GI_ratio_ohnolog_hsapiens.Pairs.Relaxed.2R.csv

generate_asymmetry_null_model_empirical_pval_ohnolog_set.R
	description: Generate null model and calculate empirical p-val comparing null model and observed ratios: ratio of negative interaction degree of paralog gene pairs from ohnolog set
	input: 
		ensembl_ohnolog_pairs_complete.rds
		Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx
		qGI_degree_auc_ess.txt
	output:
		average_null_negative_GI_ratio_seqid_50.csv
				
Run instructions:
	1. Install R and all packages listed at beginning of each script.
	2. Run 'Rscript <script_file_name>' from command line to run each script.