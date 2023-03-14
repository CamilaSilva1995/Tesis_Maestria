Los datos fueron entregados en una carpeta de drive "(https://drive.google.com/drive/folders/1x0106TYUr54gfqE6uod3g5qN8DmA_Q5x)". Para poderlos usar en el servidor, se descargaron en mi maquina local y luego se pasaron por ssh al servidor.

### En el equipo local

# $ scp Downloads/kraken_results-20230203T201634Z-001.zip camila@132.248.196.39:/home/camila/GIT/Tesis_Maestria/Data/fresa_solena
#  camila@132.248.196.39's password: 
#  kraken_results-20230203T201634Z-001.zip         100%   12MB   2.2MB/s   00:05    
# $ scp Downloads/fastp_results-20230203T175936Z-001.zip camila@132.248.196.39:/home/camila/GIT/Tesis_Maestria/Data/fresa_solena
#  camila@132.248.196.39's password: 
#  fastp_results-20230203T175936Z-001.zip          100% 2728KB   2.1MB/s   00:01  
# $ scp Downloads/bracken_results-20230203T203724Z-001.zip camila@132.248.196.39:/home/camila/GIT/Tesis_Maestria/Data/fresa_solena
#  camila@132.248.196.39's password: 
#  bracken_results-20230203T203724Z-001.zip        100% 3798KB   2.4MB/s   00:01

### En el servidor

# $ unzip kraken_results-20230203T201634Z-001.zip 
# $ unzip bracken_results-20230203T203724Z-001.zip 
# $ unzip fastp_results-20230203T175936Z-001.zip

Luego de tener las muestras Kraken y braken en el servidor con RStudio; es necesario generar el archivo **.biom**; para esto, se debe activar el ambiente de conda **metagenomis**, **kraken-biom** es un programa ampliamente utilizado para M a partir de la salida de kraken.

Con esto se obtiene una matriz de abundancia a partir de los archivos de salida de Kraken, que permite el uso de el paquete Phyloseq en R, que nos permiten analizar la diversidad y la abundancia mediante la manipulación de datos de asignación taxonómica.

# $ conda activate metagenomics

# $ kraken-biom kraken_results/* --fmt json -o fresa_kraken.biom



