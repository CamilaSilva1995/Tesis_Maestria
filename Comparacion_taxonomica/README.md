# Comparación Taxonomica

Queremos poder comparar a cualquier nivel taxonomico el resultado de la clasificación taxonomica de distintos clasificadores.


## Ejecución

```
python3 taxonomic_comparation.py -dir PATH/TO/DATA/DIRECTORY -o PATH/TO/OUTPUT/DIRECTORY
```
## Descripción

taxonomic_comparation.py lee todos los archivos con extensión .tsv del folder indicado. Estos deberían ser salidas de Kraken o Kaiju. Posteriormente une los resultados de los arhivos y hace una comparación taxonomica. La salida es un archivo .csv donde se indica read_ID, tax_ID_kraken, tax_ID_kaiju, kaiju_lineage, kraken_lineage y consensus. Donde "consensus" es el linage coincidencia de ambos programas, en caso de haber encontrado el ID en ambos clasificadores.







