# author: Arshia Z. Hassan.

R scripts:

non_essential_synthetic_lethal_analysis.R
	Descripttion: 
		1. Genrate a file with significant negative interaction pairs (non-essential library gene) and double mutant fitness scores.
		2. Generate plots: no. of unique genes and gene pairs at different Double mutant fitness score.
	input:
		Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx
		qGI_20211111_nRed.txt
		gi_FDR_20211111_nRed.txt
	output:
		SL_information.csv
		NE_SL_3_1.5.pdf
		
non_essential_synthetic_lethal_OMIM_analysis.R
	Descripttion: 
		1. Genrate a file with non-essential significant negative interaction pairs (non-essential OMIM disease library gene) and double mutant fitness scores.
		2. Generate plots: no. of unique genes and gene pairs at different Double mutant fitness score.
	input:
		morbidmap_.txt (downloaded from OMIM on 9-30-24)
		Supplmental Table 2 - Gene info_final_dec 18_2023.xlsx
		qGI_20211111_nRed.txt
		gi_FDR_20211111_nRed.txt
	output:
		SL_information.csv
		NE_OMIM_SL_3_1.5.pdf
		
		
Run instructions:
	1. Install R and all packages listed at beginning of each script.
	2. Run 'Rscript <script_file_name>' from command line to run each script.
	