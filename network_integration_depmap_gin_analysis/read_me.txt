# author: Arshia Z. Hassan.

R scripts:

pre_process.R
	Description: 
		Pre-process DepMap dependency scores.
	input:
		Achilles_gene_effect.csv
	output:
		depmap_q2_2020_nona_mean.tsv
		
rpca_onion_pipeline_qgi.R
	Description: 
		Generate RPCA-onion normalized GIN network.
	input:
		qGI_20211111_fullFF.txt		
	output:
		snf_run_rpca_7.Rdata
		
rpca_onion_pipeline.R
	Description: 
		Generate RPCA-onion normalized DepMap network.
	input:
		depmap_q2_2020_nona_mean.tsv	
	output:
		snf_run_rpca_7_5_5.Rdata
		
Integrate GIN and DepMap rpca-normalized networks with BIONIC ** 
	To generate the BIONIC integration, use settings from the GIN-DepMap-params4.json file and run code from https://github.com/bowang-lab/BIONIC.
	output:
		Gin-DepMap-params4_features.tsv

integrated_network_analysis_go_bp.r
	Description: 
		Generate PR curves and AUPRC scatter plots with GO-BP as standard.
	input:
		Mitochondial_genelist_1_26_2021_genes.tsv
		depmap_q2_2020_nona_mean.tsv
		Gin-DepMap-params4_features.tsv
		qGI_20211111_fullFF.txt
		DepMap_essential_20Q2_60_percent.txt
	output:
		AUPRC_depmap20q2_go.txt
		AUPRC_qGI_20211111_fullFF_go_bp.txt
		AUPRC_Gin-DepMap-params4_features_go_bp.txt
		AUPRC_merged_depmap_gin_bionic_rpco_go_bp.csv
		scatter_auprc_depmap_bionic_go_bp.pdf
		scatter_auprc_gin_bionic_go_bp.pdf
		scatter_auprc_depmap_gin_go_bp.pdf
		
integrated_network_analysis_corum_complex.r
	Description: 
		Generate PR curves and AUPRC scatter plots with CORUM-complex as standard.
	input:
		Mitochondial_genelist_1_26_2021_genes.tsv
		depmap_q2_2020_nona_mean.tsv
		Gin-DepMap-params4_features.tsv
		qGI_20211111_fullFF.txt
		DepMap_essential_20Q2_60_percent.txt
	output:
		AUPRC_depmap20q2_go.txt
		AUPRC_qGI_20211111_fullFF_complex.txt
		AUPRC_Gin-DepMap-params4_features_complex.txt
		AUPRC_merged_depmap_gin_bionic_rpco_complex.csv
		scatter_auprc_depmap_bionic_complex.pdf
		scatter_auprc_gin_bionic_complex.pdf
		scatter_auprc_depmap_gin_complex.pdf

Run instructions:
	1. Install R and all packages listed at beginning of each script.
	2. Run 'Rscript <script_file_name>' from command line to run each script. Run scripts in the listed order.
	3. Download FLEX package and store files in a directory to load using Devtools.
	4. Modify file paths in the scripts if necessary.

	