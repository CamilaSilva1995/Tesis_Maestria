#!/usr/bin/env python3
# coding: utf-8

# # Proyecto comparación de Outputs de clasificación taxonómica

import argparse
import pandas as pd 
import numpy as np

# lectura de los datos por medio de pandas
# Se tomaron las 10 primeras lineas de el output de kraken
datakraken10 = pd.read_csv('data/10lineas_kraken.tsv',sep='\t',header=None, names=["c-clasificado", "read_ID","tax_ID_kraken","length_pb","LCA_taxonomic"])
#datakraken10 = pd.read_csv('https://docs.google.com/spreadsheets/d/e/2PACX-1vRHYe0R7d5WYmKNbyo13egj8ORiKqWwI1GEVaG5Nzp7xN3RjkIV6vnHUweKYpidborNBf4-IEWbk0ZD/pub?output=csv',header=None, names=["c-clasificado", "read","taxon_kraken","1","2"])
# Se tomaron las 10 primeras lineas de el output de kaiju
datakaiju10 = pd.read_csv('data/10lineas_kaiju.tsv',sep='\t',header=None, names=["c-clasificado", "read_ID","tax_ID_kaiju"])
#datakaiju10 = pd.read_csv('https://docs.google.com/spreadsheets/d/e/2PACX-1vRuQ78ZkMRkMz61t0mbMBf4GGsd1ydGxXCl7p2Y8UzgXflkK6w4RCr6lPWXc9oxUXq94qeUJKouCMHp/pub?output=csv',header=None, names=["c-clasificado", "read","taxon_kaiju"])

# Se tomaron las 10 ultimas lineas de el output de kraken
datakraken10c = pd.read_csv('data/10lineascola_kraken.tsv',sep='\t',header=None, names=["c-clasificado", "read_ID","tax_ID_kraken","length_pb","LCA_taxonomic"])
#datakraken10c = pd.read_csv('https://docs.google.com/spreadsheets/d/e/2PACX-1vQJ7psycVzaws-XI_z1itUnxCI3eXBGqK4q57B4gE1YYCVCn_GkvNUZClHPMtrHANy3JUIDpqaW_r0u/pub?output=csv',header=None, names=["c-clasificado", "read","taxon_kraken2","1","2"])
# Se tomaron las 10 ultimas lineas de el output de kaiju
datakaiju10c = pd.read_csv('data/10lineascola_kaiju.tsv',sep='\t',header=None, names=["c-clasificado", "read_ID","tax_ID_kaiju2"])
#datakaiju10c = pd.read_csv('https://docs.google.com/spreadsheets/d/e/2PACX-1vQ-etvzQ8vKI0Xwggc_RvqWboa39GcRdMBw_dSa4-_mx5k9PvaMtfuaqUGvpBM3zUKLDbVkiJ3udybC/pub?output=csv',header=None, names=["c-clasificado", "read","taxon_kaiju2"])

# identificar tipo de datos**(Kraken o Kaiju) para automatizar

# Eliminamos las columnas no necesarias
datakraken101 = datakraken10.drop(["c-clasificado","length_pb","LCA_taxonomic"], axis=1)
datakaiju101 = datakaiju10.drop(["c-clasificado"], axis=1)

datakraken102 = datakraken10c.drop(["c-clasificado","length_pb","LCA_taxonomic"], axis=1)
datakaiju102 = datakaiju10c.drop(["c-clasificado"], axis=1)

# Funcion para los parametros de entrada

lOutputs = [datakraken101, datakaiju101, datakraken102, datakaiju102] # lista de dataframes

# Funcion para unir los outputs en un solo dataframe, tomando como referencia la columna 'read'
def unite(lOutputs):
    DFfusion = lOutputs[0]
    for i in lOutputs:
        DFfusion = pd.merge(DFfusion, i, how='outer')
    return DFfusion

DFfusion = unite(lOutputs)    

### Tomando la totalidad de los datos

datakraken = pd.read_csv('data/21_10_datos_completos/EG1BMAA_kraken.tsv',sep='\t',header=None, names=["c-clasificado", "read_ID","tax_ID_kraken","length_pb","LCA_taxonomic"]) 
datakaiju = pd.read_csv('data/21_10_datos_completos/EG1BMAA_kaiju.tsv',sep='\t',header=None, names=["c-clasificado","read_ID","tax_ID_kaiju"])

datakrakenall = datakraken.drop(["c-clasificado","length_pb","LCA_taxonomic"], axis=1)
datakaijuall = datakaiju.drop(["c-clasificado"], axis=1)

allOutputs = [datakrakenall, datakaijuall] #lista de dataframe
DFfusion = unite(allOutputs)    

DFfusion['measure'] = np.where((DFfusion["tax_ID_kraken"] == DFfusion["tax_ID_kaiju"]) , 1, np.nan)

DFfusion.to_csv("/data/data_ignore/clasificador.tsv", sep="\t", index=False, header=True)

        
#parser = argparse.ArgumentParser()

#parser.add_argument("-r","--read",type=str, help="read",required=True)
#parser.add_argument("-i","--pocision_inicial",type=int, help="pocision inicial", required=True)
#parser.add_argument("-n","--longitud_final",type=int, help="longitud del read", required=True)

#args = parser.parse_args()

#print()

