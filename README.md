These scripts are used for analyzing Tn-seq or RNA-seq data in the paper entitled "Widespread fungal-bacterial competition for magnesium lowers antibiotic susceptibility."

1. To analyze Tn-seq data:
   TnSeq2_copy3.sh is to map Tn-seq reads to the PAO1 genome
   run_Tnseq2_copy3.sh is to run Tnseq2_copy3.sh in BASH scripts
   Tnseq_code_for_Phoebe_2022_Nov.R is to organize the output file of TnSeq2_copy3.sh into a matrix file for the downstream DESeq2 analysis
   DESeq2_20230410_alll_coculture_test2.R is to analyze genes with depleted or enriched Tn insertions using DEseq2


2. To analyze RNA-seq data,
   20230516_RNAseq_analysis_new.R is to analyze differentially expressed genes using DEseq2   
