# author: Arshia Z. Hassan.

R scripts:

generate_all_significant_positive_GIs_suppression_metric.r
	Description: 
		1. Generate a file with significant positive GIs and suppression score.
	input:
		Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx
		qGI_20211111_nRed.txt
		gi_FDR_20211111_nRed.txt
	output:
		all_significant_positive_GIs.csv
		all_significant_positive_GIs_supp_metric_2_.9_.csv

suppressor_depmap_ess_analysis.R
	Description: 
		1. Generate boxplots of fraction of esssential in depmap for various categories of essentiality and suppression.
	input:
		Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx
		all_significant_positive_GIs_supp_metric_2_.9_.csv
	output:
		boxplot_fraction_of_esssential_depmap_suppressor_0.2.pdf
		boxplot_fraction_of_esssential_depmap_suppressor_0.5.pdf
		boxplot_fraction_of_esssential_depmap_suppressor_0.8.pdf

suppressor_depmap_ess_analysis_no_mito.R
	Description: 
		1. Generate boxplots of fraction of esssential in depmap for various categories of essentiality and suppression,
			but filter out mitochondrial genes.
	input:
		Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx
		all_significant_positive_GIs_supp_metric_2_.9_.csv
		Mitochondial_genelist_1_26_2021_genes.tsv
	output:
		boxplot_fraction_of_esssential_depmap_suppressor_no_mito_0.2.pdf
		boxplot_fraction_of_esssential_depmap_suppressor_no_mito_0.5.pdf
		boxplot_fraction_of_esssential_depmap_suppressor_no_mito_0.8.pdf

generate_all_tested_gi_data_with_suppression_score.r
	Description: 
		1. Generate a file with all tested GIs and suppression score.
	input:
		all_significant_positive_GIs_supp_metric_2_.9_.csv
		qGI_20211111_nRed.txt
		gi_FDR_20211111_nRed.txt
	output:
		all_tested_GIs_with_supp_score_metric_2_.9.csv
				
generate_PR_GO_data_flex_suppression_score.r
	Description: 
		1. Generate Precision-Recall data based on suppression score wrt to GO-BP standard using FLEX.
	input:
		Mitochondial_genelist_1_26_2021_genes.tsv
		all_tested_GIs_with_supp_score_metric_2_.9.csv
	output:
		suppressors_all_pairs_GO_metric_2_.9.Rdata
		suppressors_no_mito_GO_metric_2_.9.Rdata

generate_PR_PW_data_flex_suppression_score.r
	Description: 
		1. Generate Precision-Recall data based on suppression score wrt to pathway standard using FLEX.
	input:
		Mitochondial_genelist_1_26_2021_genes.tsv
		all_tested_GIs_with_supp_score_metric_2_.9.csv
	output:
		suppressors_all_pairs_PW_metric_2_.9.Rdata
		suppressors_no_mito_PW_metric_2_.9.Rdata

generate_PR_data_and_curve.R
	Description: 
		1. Generate PR-curve data.
	input:
		suppressors_all_pairs_GO_metric_2_.9.Rdata
		suppressors_no_mito_GO_metric_2_.9.Rdata
		suppressors_all_pairs_PW_metric_2_.9.Rdata
		suppressors_no_mito_PW_metric_2_.9.Rdata
	output:
		all_pairs_PW_PR.csv
		all_pairs_GO_PR.csv
		no_mito_PW_PR.csv
		no_mito_GO_PR.csv

generate_PR_FE_boxplot.R
	Description: 
		1. Generate barplots of AUPRC fold enrichment for different suppression score cut-offs generated against GO and Pathway standards.
	input:
		all_pairs_PW_PR.csv
		all_pairs_GO_PR.csv
		no_mito_PW_PR.csv
		no_mito_GO_PR.csv
	output:
		all_pairs_PW_PR_FE_barplot.pdf
		all_pairs_PW_PR_FE_data.csv	
		all_pairs_PW_PR_fisher_exact_results.csv
		all_pairs_GO_PR_FE_barplot.pdf
		all_pairs_GO_PR_FE_data.csv				
		all_pairs_GO_PR_fisher_exact_results.csv
		no_mito_PW_PR_FE_barplot.pdf
		no_mito_PW_PR_FE_data.csv
		no_mito_PW_PR_fisher_exact_results.csv
		no_mito_GO_PR_FE_barplot.pdf
		no_mito_GO_PR_FE_data.csv
		no_mito_GO_PR_fisher_exact_results.csv
			
Run instructions:
	1. Install R and all packages listed at beginning of each script.
	2. Run 'Rscript <script_file_name>' from command line to run each script. Run scripts in the listed order.
	3. Download FLEX package and store files in a directory to load using Devtools.
	4. Modify file paths in the scripts if necessary.
	

