import argparse
import pandas as pd
import os
from ete3 import NCBITaxa
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-dir', '--directory', help='data path')
parser.add_argument('-o', '--output', help = 'output path')
args = parser.parse_args()


PATH = args.directory


def read_files(directory):
    kraken, kaiju = [], []
    for file in os.listdir(directory):
        # print(file)
        if '.tsv' in file:
            f = pd.read_csv(directory + file, delimiter = '\t', header = None)
            if f.shape[1] == 3:
                f.columns = ["c-clasificado","read_ID","tax_ID_kaiju"]
                f.drop(['c-clasificado'], axis=1, inplace = True)
                kaiju.append(f)

            else:
                f.columns = ["c-clasificado", "read_ID","tax_ID_kraken","length_pb","LCA_taxonomic"]
                f.drop(["c-clasificado","length_pb","LCA_taxonomic"], axis=1, inplace = True)
                kraken.append(f)
            
    kraken, kaiju = pd.concat(kraken), pd.concat(kaiju)
        
    return [kraken, kaiju]

def unite(lOutputs):
    DFfusion = lOutputs[0]
    for i in lOutputs:
        DFfusion = pd.merge(DFfusion, i, how='outer')
   
    return DFfusion


def get_full_lineages(otus):

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
        lineage = ncbi.get_lineage(entrie)               
        names = ncbi.get_taxid_translator(lineage)
        names_ordered = []
        for l in lineage:
            names_ordered.append(names[l])
        
        all_names = ";".join(names_ordered)
        lineages.update({entrie: all_names})
    
    return lineages
    

def add_lineages(DFfusion):
    kaiju_lineages = get_full_lineages(list(DFfusion['tax_ID_kaiju'].dropna()))
    kraken_lineages = get_full_lineages(list(DFfusion['tax_ID_kraken'].dropna()))
    DFfusion['kaiju_lineage'] = DFfusion['tax_ID_kaiju']
    DFfusion['kraken_lineage'] = DFfusion['tax_ID_kraken']
    DFfusion.replace({'kaiju_lineage':kaiju_lineages, 'kraken_lineage':kraken_lineages}, inplace = True)
    
    return DFfusion

def consensus(lineages_df):
    l = []
    for item in zip(lineages_df['kraken_lineage'].fillna(''), lineages_df['kaiju_lineage'].fillna('')):
        if item[0] != '' and item[1] != '':
            lineage_a = item[0].split(';')
            lineage_b = item[1].split(';')
            lineage_intersection = []
            if len(lineage_a) < len(lineage_b):
                for i in range(len(lineage_a)):
                    if lineage_a[i] == lineage_b[i]:
                        lineage_intersection.append(lineage_a[i])
                    else:
                        break
            else:
                for i in range(len(lineage_b)):
                    if lineage_a[i] == lineage_b[i]:
                        lineage_intersection.append(lineage_b[i])
                    else:
                        break

            l.append(';'.join(lineage_intersection))
        else:
            l.append(np.nan)
        
    lineages_df['consensus'] = l
    return lineages_df



outputs = consensus(add_lineages(unite(read_files(PATH))))

if not os.path.isdir(args.output):
    os.makedirs(args.output)

if args.output[-1] == '/':
    outputs.to_csv(args.output + 'results.csv')
else:
    outputs.to_csv(args.output + '/results.csv')



