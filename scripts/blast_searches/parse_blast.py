#Parse blast results and return top hit tsv and fasta file of sequences

import csv
import sys
import Bio
import os
import re

from Bio import SeqIO

def search_fasta(input_fasta, output_fasta, results): 
    seq_records = SeqIO.parse(input_fasta, format='fasta') #parses the fasta file
    
    with open(results, 'r') as file:
        top_hits = csv.reader(file, delimiter='\t')
        next(top_hits)
        ids_to_get = []
        for hit in top_hits:
            ids_to_get.append(hit[1])
    
    with open(output_fasta, 'w') as outfile:
        for record in seq_records: 
            if record.id in ids_to_get: 
                outfile.write(f">{record.id}\n")
                outfile.write(f"{record.seq}\n")
			


def search_fasta_un(input_fasta, output_fasta, results): 
    seq_records = SeqIO.parse(input_fasta, format='fasta') #parses the fasta file
    
    with open(results, 'r') as file:
        top_hits = csv.reader(file, delimiter='\t')
        next(top_hits)
        ids_to_get = []
        coords = {}
        for hit in top_hits:
            ids_to_get.append(hit[1])
            coords[hit[1]] = [hit[-4], hit[-3]]
    
    with open(output_fasta, 'w') as outfile:
        for record in seq_records: 
            if record.id in ids_to_get: 
                if int(coords[record.id][0]) <= int(coords[record.id][1]):
                    start = int(coords[record.id][0]) 
                    end = int(coords[record.id][1]) 
                else:
                    end = int(coords[record.id][0]) 
                    start = int(coords[record.id][1]) 
                outfile.write(f">{record.id}_{start}_{end}\n")
                outfile.write(f"{record.seq[start:end+1]}\n")
			



def parse_blast_results(input, output):

    with open(input) as file:
        contents = csv.reader(file, delimiter='\t')
        top_hits = {}
        for gene in contents:
            if gene[0] not in top_hits.keys():
                top_hits[gene[0]] = gene

    print(f"Number of hits: {len(top_hits)}") 

    with open(output, 'w') as result:
        result.write(f"gene\thit\tpercentid\taln_length\tsub-start\tsub-end\tevalue\tbitscore\n")
        for key,value in top_hits.items():
            result.write(f"{key}\t{value[1]}\t{value[2]}\t{value[3]}\t{value[-4]}\t{value[-3]}\t{value[-2]}\t{value[-1]}\n")

    
        

def main():

    try:
        os.mkdir("tophits_results")
        os.mkdir("tophits_fasta")
    except OSError as error:
        print(error)    
    
    
    for file in os.listdir("results"):
        basename = re.split("_results", file)[0]
        print(f"{basename}")
        parse_blast_results(f"results/{file}", f"tophits_results/{basename}.tophits.tsv")
        if os.path.isfile(os.path.join('annotated', f"{basename}_protein.faa")):
            search_fasta(f"annotated/{basename}_protein.faa", f"tophits_fasta/{basename}.tophits.faa", f"tophits_results/{basename}.tophits.tsv")
        else:
            search_fasta_un(f"unannotated/{basename}.fna", f"tophits_fasta/{basename}.tophits.fna", f"tophits_results/{basename}.tophits.tsv")

if __name__ == "__main__":
    main()

