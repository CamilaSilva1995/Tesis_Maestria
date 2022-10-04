#!/usr/bin/env python3
# coding: utf-8

import argparse, os               
import pandas as pd  

### Funcion para recortar reads

# Funcion cutout(string a recortar, pocision de inicio, longitud del recorte)
def cutout(read,i,n_length): 
    cropped_read = read[i:i+n_length]
    return(cropped_read)
    
    #queremos un Output como archivo .tsv que contenga (pocision_inicial, longitud_final, reads_cortos)
    #reads_cortos = pd.DataFrame(cropped_read.items(), columns=["i", "n", "read"])
    #reads_cortos.to_csv("reads_cortos.tsv", sep="\t", index=False, header=True)

# ejemplo  
#cropped_read = cutout("CCGTTCCTCGGGCGTGCAGTCGGGCTTGCGGTCTGCCATGTCGTGTTCGGCGTCGGTGGTGCCGATCAGGGTGAAATCCGTCTCGTAGGGGATCGCGAAGATGATCCGCCCGTCCGTGCCCTGAAAGAAATAGCACTTGTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCTCAGAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAGCAAACCTCTCACTCCCTCTACTCTACTCCCTT",25,10)
#cropped_read

####### Arguments #######
    
parser = argparse.ArgumentParser()

parser.add_argument("-r","--input_archivoFasta", help="read",required=True)
parser.add_argument("-i","--pocision_inicial", type=int, help="pocision inicial", required=True)
parser.add_argument("-n","--longitud_final", type=int, help="longitud del read", required=True)

args = parser.parse_args()

print(cutout(args.read, args.pocision_inicial, args.longitud_final))



