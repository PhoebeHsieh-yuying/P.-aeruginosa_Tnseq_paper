#!/bin/bash
#
#SBATCH --job-name run_Tnseq2_copy3
#SBATCH --cpus-per-task=4
#SBATCH --time=120:00:00
#SBATCH --partition=campus-new
#SBATCH --mem-per-cpu=5g

module load cutadapt
#module load Python
module load Bowtie2

gzip -d Mg_1_CKDL220032871-1A_HJCYKBBXX_L3_1.fq.gz

mv Mg_1_CKDL220032871-1A_HJCYKBBXX_L3_1.fq
Mg_1_CKDL220032871-1A_HJCYKBBXX_L3_1.fastq

/fh/fast/malik_h/user/yhsieh/Raw_Sequence_Reads/Tnseq_20221226/TnSeq2_copy3.sh -i TATAAGAGTCAG -g /fh/fast/malik_h/user/yhsieh/PAO1_Reference_Genome/Bowtie2_indice/PAO1 Mg_1_CKDL220032871-1A_HJCYKBBXX_L3_1
