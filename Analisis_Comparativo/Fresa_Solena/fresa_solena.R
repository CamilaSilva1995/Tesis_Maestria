### DIVERSIDADES CON TODO EL CONJUNTO DE DATOS  Y LOS DATOS FILTRADOS POR CALIDAD
## Los datos fueron entregados en una carpeta de drive...
## para poderlos usar en el servidor, se descargaron en mi maquina local y luego los pase por ssh al servido 

## Local
# $ scp Downloads/kraken_results-20230203T201634Z-001.zip camila@132.248.196.39:/home/camila/GIT/Tesis_Maestria/Data2/fresa_solena
# $ camila@132.248.196.39's password: 
# $ kraken_results-20230203T201634Z-001.zip         100%   12MB   2.2MB/s   00:05    

# $ scp Downloads/fastp_results-20230203T175936Z-001.zip camila@132.248.196.39:/home/camila/GIT/Tesis_Maestria/Data2/fresa_solena
# $ camila@132.248.196.39's password: 
# $ fastp_results-20230203T175936Z-001.zip          100% 2728KB   2.1MB/s   00:01    

# $ scp Downloads/bracken_results-20230203T203724Z-001.zip camila@132.248.196.39:/home/camila/GIT/Tesis_Maestria/Data2/fresa_solena
# $ camila@132.248.196.39's password: 
# $ bracken_results-20230203T203724Z-001.zip        100% 3798KB   2.4MB/s   00:01

## En el servidor
# $ unzip kraken_results-20230203T201634Z-001.zip 
# $ unzip bracken_results-20230203T203724Z-001.zip 
# $ unzip fastp_results-20230203T175936Z-001.zip

# $ conda activate metagenomics
# $ kraken-biom kraken_results/* --fmt json -o fresa_kraken.biom

##
getwd()
## Es necesario instalar phyloseq de la siguiente manera 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  +     install.packages("BiocManager")
#BiocManager::install("phyloseq") 

#install.packages(c("ggplot2", "readr", "patchwork"))
## Ahora liego de instalar las librerias, necesitamos llamarlas a la sesion de R en la que se va a trabajar.
library("phyloseq") 
## Phyloseq es un paquete de Bioconductor (Open Source Software For Bioinformatics) para la manipulación y análisis (herramienta para importar, guardar, analizar y visualizar) de datos metagenómicos generados por metodologías de secuenciación de alto rendimiento. 
library("ggplot2") #
library("igraph")
library("readr")
library("patchwork")
library("vegan")
library("GUniFrac")
library("pbkrtest")
#library("BiodiversityR")
library("kableExtra")
library("RColorBrewer")
## Debemos especificarle a R en que directorio estamos trabajando
setwd("/home/camila/GIT/Tesis_Maestria/Data2/fresa_solena")

## Datos kraken
## importamos los datos kraken del archivo biom, usando la herramienta de phyloseq **impor_biom**
fresa_kraken <- import_biom("fresa_kraken.biom")
class(fresa_kraken) # objeto phyloseq

## Queremos acceder a los datos que contiene nuestro objeto phyloseq **fresa_kraken**
## primero la tabla de taxonomia
#TAX = fresa_kraken@tax_table@.Data
## cambiamos los nombres de las columnas a los niveles taxonomicos
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
## Queremos cortar la parte inicial nombres, ya que aparecen,por ejemplo: "B__Bacteria" y queremos que solo se vea "Bacteria"
## substring(). Este comando ayuda a extraer o reemplazar caracteres en un vector.
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data,4,100)
#View(fresa_kraken@tax_table@.Data)
## Queremos ver la tabla de OTUs, esta tabla contiene las abundancias de los otus de cada muestra
#OTU = fresa_kraken@otu_table@.Data
## Queremos cortar los nombres de las muestras para que coincida con los metadatos 
colnames(fresa_kraken@otu_table@.Data) <- substr(colnames(fresa_kraken@otu_table@.Data),1,6)
#View(fresa_kraken@otu_table@.Data)

### cargar los metadatos 
### ya que hay un desface de dos muestras entre los metadatos y las muestras de la otu_table
### se deben ver cuales son y quitarlas del archivo de metadatos

# $ ls kraken_results |cut -d'.' -f1 > lista_kraken.txt
# $ ls metadata.csv |cut -d',' -f1 > lista_metadata.txt
# $ wc *txt
# $ cat lista_metadata.txt lista_kraken.txt | sort | uniq -c
# $ cat lista_metadata.txt lista_kraken.txt | sort | uniq -c | sort | head
# $  1 MD2145
# $  1 MD2146
# $  2 MD2055
# $  2 MD2056

## ELIMINAR LOS DATOS QUE SOBRAN DE LOS METADATOS
### LO HICEA MANO... REVISAR COMO HACERLO EN BASH
## eliminamos las dos muestras que no estaban en nuestra otu_table y cargamos los metadatos
metadata_fresa <- read.csv2("/home/camila/GIT/Tesis_Maestria/Data2/fresa_solena/metadata.csv",header =  FALSE, row.names = 1, sep = ",")
#metadata_fresa <- sample_data(metadata_fresa)
## luego hacemos que los metadatos pertenezcan al objeto phyloseq en la seccion de **sam_data**
fresa_kraken@sam_data <- sample_data(metadata_fresa) 
## Creamos una columna extra en sam_data por necesidad de funcionamiento de mas adelante
fresa_kraken@sam_data$Sample<-row.names(fresa_kraken@sam_data)
colnames(fresa_kraken@sam_data)<-c('Treatment','Samples')
#View(fresa_kraken)

## Ahora tenemos tambien la tabla de datos de calidad (fastp_kraken_summary), 
## para ver que muestra podemos eliminar de nuestro dataset, que no cumpla ciertos estandares de calidad
## ID de la muestra 
## Reads_B - Reads_Before -> total de reads crudos
## Reads_A - Reads_After -> total de reads despues del analisis de calidad
## Reads_diff -> diferencia en tre Reads_B y Reads_A 
## Q30_B -> porcentaje arriba de 30 (escala fred) antes del analisis de calidad 
## Q30_A -> porcentaje arriba de 30 (escala fred) despues del analisis de calidad
## LowQua -> reads de baja calidad
## N_reads -> readas que contienen N y se descartan
## too_short -> no pasan el tamaño minimo de calidad
## Duplication -> porcentaje de duplicados
## LengthR1 -> longitud promedio de los reads 
## LengthR2 -> longitud promedio de los reads 
## Classified -> porcentaje de clasificados del total despues del filtrado 

## ejemplo muestra MD2055 -> contienen 97millones de reads antes del filtrado de calidad, y despues queda con 79millones
## filtro usado para eliminar muestras, que luego del filtrado de calidad contengan menos de 25millones de reads 
## los que nos da 5 muestras a eliminar (MP2079,MP2080,MP2088,MP2109,MP2137)

## eliminiar muestras de baja calidad,usando filtro de menos de 25 millones de reads luego del analisis de calidad
## se eliminan las muestras por su nombre
samples_to_remove <- c("MP2079","MP2080","MP2088","MP2109","MP2137")
fresa_kraken_fil <- prune_samples(!(sample_names(fresa_kraken) %in% samples_to_remove), fresa_kraken)
## podemos comprobar el número de muestras antes y después del filtrado
nsamples(fresa_kraken) # 58
nsamples(fresa_kraken_fil) # 53

## Queremos hacer un analisis de diversidad de nuestras muestras, para esto las dos metricas las usadas son: Diversidad Alfa y Beta

## Diversidad Alfa
## Esta representa la riqueza de las muestras, es decir el numero de especies diferentes en ese ambiente o la abundancia de especies en ese ambiente. Paramedier esta divercidad se tienen diferentes indices de medida, entre ellos losindices Shannon,Simpson, Chao1, entre otros.
index = estimate_richness(fresa_kraken_fil)
## podemos ver una representación visual de la diversidad dentro de las muestras
## esta diversidad la podemos ver por muestra,
plot_richness(physeq = fresa_kraken_fil, measures = c("Observed","Chao1","Shannon","simpson"))
ggsave("DiversidadAlfa.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
## pero ya que queremos diferenciar entre plantas sanas y enfermas, tomamos estos dos conjuntos por ceparado para mostrar la diversidad
plot_richness(physeq = fresa_kraken, measures = c("Observed","Chao1","Shannon","simpson"),x = "Treatment", color = "Treatment") 
ggsave("DiversidadAlfaFresaKraken.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img", width = 30, height = 15, dpi = 300, units = "cm")
## comparando entre los datos en crudo,y los datos filtrados por calidad
plot_richness(physeq = fresa_kraken_fil, measures = c("Observed","Chao1","Shannon","simpson"),x = "Treatment", color = "Treatment") 
ggsave("DiversidadAlfaFresaKraken_fil.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

## Shannon
##
## Simpson
##
## Chao1
##

###############################################################################################################################

## Aqui podemos ver cuantos reads tenemos por muestra
sample_sums(fresa_kraken)
sample_sums(fresa_kraken_fil)
## y con **summary** podemos darnos una idea dela uniformidad de los datos, ya que podemos ver datos estadisticos como elmaximo, el minimo y la media
summary(fresa_kraken@otu_table@.Data)
summary(fresa_kraken_fil@otu_table@.Data)

## Con el siguiente comando podemos ver si tenemos muestras no identificadas taxonomicamente, esto se puede ver identificando los espacios en blanco ("") en los diferentes niveles taxonomicos.
summary(fresa_kraken_fil@tax_table@.Data== "") 
## podemos ver los **TRUE** de cada nivel taxonomico,por ejemplo a nivel de "Phylum" tenemos solo 2 sin clasificar, y a nivel de "Specie" tenemos 1099 sin clasificar
## ESTO ES LO QUE QUEREMOS COMPARAR CON LAS SALIDAS DE BRACKEN, YA QUE PROMETE HACER UNA REASISGNACION DE TODO LO QUE QUEDE SIN CLASIFICAR CON KRAKEN

## Queremos convertir las abundancias absolutas a relativas, Calculamos las abundancias relativas con la función x/sum(x)∗100, tanto de los datos originales como delos filtrados 
percentages <- transform_sample_counts(fresa_kraken, function(x) x*100 / sum(x) )
head(percentages@otu_table@.Data)
percentages_fil <- transform_sample_counts(fresa_kraken_fil, function(x) x*100 / sum(x) )
head(percentages_fil@otu_table@.Data)

## Ahora ya es posibleuna buena comparacion de abundancias dadas por porcentajes con indices beta 

## Diversidad Beta
## La diversidad beta mide la diferencia entre dos o mas entornos. Se puede medir con métricas como la disimilitud de Bray-Curtis, la distancia Jaccard o la distancia UniFrac.
## Mide que tan similares o diferentes son un par de especies, muestras o conjuntos de muestras.
## aqui podemos ver una lista de distancias disponibles, que Phyloseq puede usar
distanceMethodList 
## Siendo las siguientes las mas usadas:
## Disimilitud de Bray-Curtis
##
## Distancia Jaccard
##
## UniFrac
##

vegdist(meta_ord_fil,"bray")# usa lalibreria Vegan,  REVISAR COMO FUNCIONA PARA IMPRIMIR LA TABLA

## Usamos "ordinate" para asignar las distancias entre muestras,usando "Bray-Curtis", ya que es una de las metricas mas completas y mayormente utilizadas para medir la diversidad beta
## Hay diferentes formas de trazar y mostrar los resultados de dicho análisis. Entre otros, se utilizan ampliamente los análisis PCA, PCoA o NMDS.
meta_ord <- ordinate(physeq = percentages, method = "NMDS", distance = "bray")   
meta_ord_fil <- ordinate(physeq = percentages_fil, method = "NMDS", distance = "bray") 
## Ademas lo queremos diferenciar por color entre plantas sanas y enfermas,

#####ggsave("DiversidadAlfaFresaKraken_fil.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

plot_ordination(physeq = percentages, ordination = meta_ord, color = "Treatment") + 
  geom_text(mapping = aes(label = colnames(fresa_kraken@otu_table@.Data)), size = 3, vjust = 1.5)
ggsave("DiversidadBetaFresaKraken.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
plot_ordination(physeq = percentages_fil, ordination = meta_ord_fil, color = "Treatment") +
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)
ggsave("DiversidadBetaFresaKraken_fil.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
## En este gráfico NMDS, cada punto representa la abundancia combinada de todas sus OTU. No se puiede ver una diferencia entre los conjuntos de muestras de plantas sanas y muestras enfermas, por lo tanto:
## Probaremos varias distancias ya con los datos filtrados podemos ver varios ejemplos y llegar  a la posibilidad que se diferencien un poco los dos conjuntos de datos
## Usando 
meta_ord_fil_2 <- ordinate(physeq = percentages_fil, method = "NMDS", distance = "jsd")
plot_ordination(physeq = percentages_fil, ordination = meta_ord_fil_2, color = "Treatment") + 
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)
ggsave("DiversidadBetaFresaKraken_fil_jsd.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
## Usando 
meta_ord_fil_3 <- ordinate(physeq = percentages_fil, method = "NMDS", distance = "euclidean")
plot_ordination(physeq = percentages_fil, ordination = meta_ord_fil_3, color = "Treatment") + 
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)
ggsave("DiversidadBetaFresaKraken_fil_euclidean.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
## Usando 
meta_ord_fil_4 <- ordinate(physeq = percentages_fil, method = "NMDS", distance = "jaccard")
plot_ordination(physeq = percentages_fil, ordination = meta_ord_fil_4, color = "Treatment") + 
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)
ggsave("DiversidadBetaFresaKraken_fil_jaccard.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
## Usando 
meta_ord_fil_5 <- ordinate(physeq = percentages_fil, method = "NMDS", distance = "manhattan")
plot_ordination(physeq = percentages_fil, ordination = meta_ord_fil_5, color = "Treatment") + 
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)
ggsave("DiversidadBetaFresaKraken_fil_manhattan.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
## Con ninguna distancia vemos una diferenciacion clara, por lo que lo veremos por los distintos niveles taxonomicos 



