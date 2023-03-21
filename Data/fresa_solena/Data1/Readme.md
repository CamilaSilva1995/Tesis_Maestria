Los datos fueron entregados en una carpeta de drive https://drive.google.com/drive/folders/1x0106TYUr54gfqE6uod3g5qN8DmA_Q5x. Para poderlos usar en el servidor, se descargaron en mi maquina local y luego se pasaron por ssh al servidor.

### En el equipo local
```{bash}
$ scp Downloads/kraken_results-20230203T201634Z-001.zip camila@132.248.196.39:/home/camila/GIT/Tesis_Maestria/Data/fresa_solena
camila@132.248.196.39's password: 
kraken_results-20230203T201634Z-001.zip         100%   12MB   2.2MB/s   00:05    
$ scp Downloads/fastp_results-20230203T175936Z-001.zip camila@132.248.196.39:/home/camila/GIT/Tesis_Maestria/Data/fresa_solena
camila@132.248.196.39's password: 
fastp_results-20230203T175936Z-001.zip          100% 2728KB   2.1MB/s   00:01  
$ scp Downloads/bracken_results-20230203T203724Z-001.zip camila@132.248.196.39:/home/camila/GIT/Tesis_Maestria/Data/fresa_solena
camila@132.248.196.39's password: 
bracken_results-20230203T203724Z-001.zip        100% 3798KB   2.4MB/s   00:01
```

### En el servidor

```{bash}
$ unzip kraken_results-20230203T201634Z-001.zip 
$ unzip bracken_results-20230203T203724Z-001.zip 
$ unzip fastp_results-20230203T175936Z-001.zip
```

Luego de tener las muestras Kraken y braken en el servidor con RStudio; es necesario generar el archivo **.biom**; para esto, se debe activar el ambiente de conda **metagenomis**, **kraken-biom** es un programa ampliamente utilizado para M a partir de la salida de kraken.

Con esto se obtiene una matriz de abundancia a partir de los archivos de salida de Kraken, que permite el uso de el paquete Phyloseq en R, que nos permiten analizar la diversidad y la abundancia mediante la manipulación de datos de asignación taxonómica.

```{bash}
$ conda activate metagenomics
$ kraken-biom kraken_results/* --fmt json -o fresa_kraken.biom
```


## Modificando los datos para poder usar  las herramientas de redes:

### Alnitak

132.248.196.39:5000

Calcula los vecinos de un taxon, (taxón de interes con dos metricas)
```{r}
# tomamos un subconjunto cortado a nivel de genero, objeto phyloseq
Genus_glom <- tax_glom(fresa_kraken_fil, taxrank = 'Genus')
# se pasa a dataframe
Genus_glom_df = Genus_glom@otu_table
# se guarda como csv
write.csv(Genus_glom_df, "/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/fresa_genus.csv")
```
Se tuvo que corre un lugar de losnombre delascolumnas,ya quegeneraba un error al tratar de convertir ese archivo tsv a .BIOM
```{bash}
$ biom convert -i fresa_genus.tsv -o freas_genus.biom --table-type="OTU table" --to-json
```
Se debe subir:
1. TaxID
2. Documento .BIOM cortado por nivel taxonomico del taxon escojido
Entrega matriz de correlaciones filtrada porel taxon de interes

Output:
* 20230314-1_raw_network.csv

Que sólo reporta las correlaciones con el género que se le pone y de esas, sólo si pasan los umbrales de:

* Corrección > 0.7 y
* Bray Curtis< 0.3

Que cumpla estas condiciones tomando a nivel de genero con taxon principal Fusarium (5506) tenemos: taxon2 (5125,57161,37994,5141)


### MicNet

132.248.196.39:8501

Conda MicNet -env

Se requiere:
Documento BIOM en formato csv, esto es posible usando la terminal de la siguiente manerfa 
```{bash}
$ biom convert -i fresa_kraken.biom -o fresa_kraken.tsv --to-tsv
$ sed 's/\t/,/g' fresa_kraken.tsv > fresa_kraken.csv
```
1. UMAP -> Cluster
2. Sparcc -> Matrizde correlacion
3. Network ->Pidelos documentos generados en los dos anteriores pasos


MicNet crea la red de coocurrencia completa, entregando lamatrizde correlaciones.
Usa normalizacion de Dirichlet

Output:
Se obtienen 

* Output_UMAP_HDBSCAN.csv
* SparCC_Output.csv
* 
