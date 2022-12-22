#!/usr/bin/env python3
# coding: utf-8

import argparse, os               
import pandas as pd  

def get_full_lineages(otus):

    #### makes the full lineage file (lineages.tsv). Requires ete3 ####

    #Input: list of the otus in the table obtained from the get_otus function
    #Output: makes the full_lineages.tsv file 

    from ete3 import NCBITaxa

    ncbi = NCBITaxa()  

    lineages = {}  

    if 0 in otus:
        lineages.update({0: ""})
        otus.remove(0)
    if 1 in otus:
        lineages.update({1: "root"})
        otus.remove(1)
    if 2 in otus:
        lineages.update({2: "root;Bacteria"})
        otus.remove(2)

    for entrie in otus:
        lineage = ncbi.get_lineage(entrie)                  #returns list of lineage taxids
        names = ncbi.get_taxid_translator(lineage).values() #returns dict in which the taxids of the lineage list become the keys (int) and the translations the values. Error if there is a 0 
        all_names = ";".join(names)
        lineages.update({entrie: all_names})

    #queremos que nos arroje un archivo .tsv => (READ/LINAGE1/LINAJE2, etc...)    
    lineages_df = pd.DataFrame(lineages.items(), columns=["OTU", "LINEAGE"])
    lineages_df.to_csv("full_lineages.tsv", sep="\t", index=False, header=True)
    print("full lineage file created")
    
# example
#get_full_lineages([0,329,547])

####### Arguments #######
    
parser = argparse.ArgumentParser()

parser.add_argument("-i","--input_file", help="file with otus",required=True)

args = parser.parse_args()

print(get_full_lineages(args.otus))
