#Parse blast results and return top hit tsv and fasta file of sequences

import csv
import sys
import Bio
import os
import re

from Bio import SeqIO

def parse_blast_results(input, output):

    #parse results for blast search and return top match 
    #which is always first instance in file
    #this assumes only want one top match 
    with open(input) as file:
        contents = csv.reader(file, delimiter='\t')
        top_hits = {}
        for gene in contents:
            if gene[0] not in top_hits.keys():
                top_hits[gene[0]] = gene

    print(f"Number of reciprical hits: {len(top_hits)}") 

    with open(output, 'w') as result:
        result.write(f"gene\thit\tpercentid\taln_length\tsub-start\tsub-end\tevalue\tbitscore\n")
        for key,value in top_hits.items():
            result.write(f"{key}\t{value[1]}\t{value[2]}\t{value[3]}\t{value[-4]}\t{value[-3]}\t{value[-2]}\t{value[-1]}\n")

    
def compare_blast_results(first, last, result_file):
    #compare initial blast search parsed results to final parsed results

    keyfile = get_symbiosis_trans("symbiosis_genes.txt")

    with open(first) as file1, open(last) as file2:
        initial = csv.reader(file1, delimiter='\t')
        last = csv.reader(file2, delimiter='\t')
        next(initial)
        next(last)

        matches = {}
        results = {}

        for hit in initial:
            matches[hit[1]] = [hit[0], keyfile[hit[0]]]
        
        for hit in last: 
            try:
                geneID = re.split(r'_\d*_\d*', hit[0])[0]
                matchID = re.split(r'\.\d', hit[1])[0]
                print(matches[geneID][1])
                print(matchID)
                if geneID in matches.keys():
                    if matches[geneID][1] == matchID:
                        results[geneID] = True
                    else:
                        results[geneID] = False
                else:
                    results[geneID] = 'NA'
            except:
                print('Error')       

    with open(result_file, 'w') as res:
        
        res.write(f"Gene\tHit\tRecipricalMatch\n")
        for key, val in results.items():
            res.write(f"{matches[key][0]}\t{key}\t{val}\n")



def get_symbiosis_trans(keyfile):

    with open(keyfile) as k:
        contents = csv.reader(k, delimiter='\t')
        next(contents)

        symbiosis = {}
        for gene in contents:
            symbiosis[gene[0]] = gene[1]
        
    
    return symbiosis


def main():

    try:
        os.mkdir("match_results")
        os.mkdir("recip_top_results")
    except OSError as error:
        print(error)    
    
    
    for file in os.listdir("recip_results"):
        basename = re.split("_results", file)[0]
        print(f"{basename}")
        parse_blast_results(f"recip_results/{file}", f"recip_top_results/{basename}.tophits.tsv")
        compare_blast_results(f"tophits_results/{basename}.tophits.tsv", f"recip_top_results/{basename}.tophits.tsv", f"match_results/{basename}.results.tsv")

if __name__ == "__main__":
    main()

