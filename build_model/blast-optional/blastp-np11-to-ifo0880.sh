#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --mem=1GB 
#SBATCH --time=10:00:00 
#SBATCH --partition=open
module load ncbi-blast
blastp -query uniprotkb_proteome_UP000016926_2023_11_01.fasta -subject uniprotkb_proteome_UP000239560_2023_11_01.fasta -out blastp-np11-to-ifo0880.txt -evalue 1e-5 -outfmt '6 qseqid sseqid score bitscore evalue pident qcovs sallseqid'
