#!/bin/bash -l
#
#SBATCH --ntasks 16 #number cores
#SBATCH -o logs/02_recipr.log
#SBATCH -e logs/02_recipr.log
#SBATCH -J recip_search
#SBATCH --mem=50G #memory
#SBATCH -p short


module load ncbi-blast


# make reference DB
# makeblastdb -in reference/Mtruncatula_285_Mt4.0v1.protein.fa -dbtype 'prot' -out reference/Mtruncatula_285_Mt4.0v1

CPU=16
INPUT=tophits_fasta
REFERENCE=reference/Mtruncatula_285_Mt4.0v1
RES=recip_results

mkdir $RES

for file in $(ls $INPUT/*faa);
do 
	base=$(basename $file .tophits.faa)
	echo $base
	blastp -db $REFERENCE -num_threads $CPU -outfmt 6 -evalue 1e-10 -query $file -out $RES/$base"_results"

done 

for file in $(ls $INPUT/*.fna);
do 
	base=$(basename $file .tophits.fna)
	echo $base
	blastx -db $REFERENCE -num_threads $CPU -outfmt 6 -evalue 1e-10 -query $file -out $RES/$base'_results'

done 


