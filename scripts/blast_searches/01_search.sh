#!/bin/bash -l
#
#SBATCH --ntasks 16 #number cores
#SBATCH -o logs/01_search.log
#SBATCH -e logs/01_search.log
#SBATCH -J initial_search
#SBATCH --mem=50G #memory
#SBATCH -p intel
#SBATCH --time=24:00:00


module load ncbi-blast


# make reference DB
# makeblastdb -in reference/Mtruncatula_285_Mt4.0v1.protein.fa -dbtype 'prot' -out reference/Mtruncatula_285_Mt4.0v1

CPU=16
ANNOTATED=annotated
UNANNOTATED=unannotated
REFERENCE=reference/Mtruncatula_285_Mt4.0v1
SYMBIOSIS=symbiosis_genes.fa
RES=results

for file in $(ls $ANNOTATED/*faa);
do 
	base=$(basename $file _protein.faa)
	echo $base
	makeblastdb -in $file -dbtype 'prot' -out $ANNOTATED/$base
	blastp -db $ANNOTATED/$base -num_threads $CPU -outfmt 6 -evalue 1e-10 -query $SYMBIOSIS -out $RES/$base"_results"

done 

for file in $(ls $UNANNOTATED/*.fna);
do 
	base=$(basename $file .fna)
	echo $base
	makeblastdb -in $file -dbtype 'nucl' -out $UNANNOTATED/$base
	tblastn -db $UNANNOTATED/$base -num_threads $CPU -outfmt 6 -evalue 1e-10 -query $SYMBIOSIS -out $RES/$base'_results'

done 




