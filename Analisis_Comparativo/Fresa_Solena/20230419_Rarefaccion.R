## RAREFACCION (RICHNESS)
# La curva de rarefacción da la riqueza de especies anotadas es un gráfico, el
# número total de anotaciones de especies distintas en función del número de
# secuencias muestreadas.

library("phyloseq")
library("ggplot2")
library("vegan")
library("RColorBrewer")
library("stringi")
library("dplyr")
library("plyr")

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
### Cargado de datos originales 
fresa_kraken <- import_biom("fresa_kraken.biom")
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data,4,100)
colnames(fresa_kraken@otu_table@.Data) <- substr(colnames(fresa_kraken@otu_table@.Data),1,6)
metadata_fresa <- read.csv2("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1/metadata.csv",header =  FALSE, row.names = 1, sep = ",")
fresa_kraken@sam_data <- sample_data(metadata_fresa)
fresa_kraken@sam_data$Sample<-row.names(fresa_kraken@sam_data)
colnames(fresa_kraken@sam_data)<-c('Treatment','Samples')
samples_to_remove <- c("MP2079","MP2080","MP2088","MP2109","MP2137")
fresa_kraken_fil <- prune_samples(!(sample_names(fresa_kraken) %in% samples_to_remove), fresa_kraken)
fresa_kraken_fil_df <- psmelt(fresa_kraken_fil)

## solo tomamos las columnas Sample, OTU y Abundance
fresa_kraken_fil_df2 <- fresa_kraken_fil_df[c(1,2,3)]

# queremos pasar de dataframe a table, para calcular el alfa diversidad
dfss <- reshape(fresa_kraken_fil_df2, idvar = "Sample",v.names= c("Abundance"), timevar = "OTU", direction = "wide")
rownames(dfss) <- dfss$Sample
dfss <- select(dfss, -Sample)

col<-specaccum(dfss, method = "collector")
plot(col, xlab="Muestras", ylab="Abundancia", col="blue",ylim = c(8500,9000))
points(col$richness, pch=19, col="darkblue")

## Rarefaccion basada en individuos
#Sumamos la abundancia de cada especie
N <- colSums(dfss)
#Hacemos un vector con los tamaños de muestra sobre los cuales haremos 
#la interpolación. El dato final de este vector es el tamaño total de la 
#muestra (sum(N))
subs3 <- c(seq(500, 4000, by = 500), sum(N)) 
#Ejecutamos la rarefacción
rar3 <- rarefy(N, sample = subs3, se = T, MARG = 2)
rar3
# muestra la cantidad deespecies que se espera tener para diferentes tamaños de muestras

## Rarefaccion basada en muestras
#acomulacion de especies
rand <- specaccum(dfss, method = "random", permutations=100)
rand
#vemos la acomulacion de otus encontrados dependiendo la cantidad de muestras


### Rarefy (sanos y enfermos)

Treatment <-split(fresa_kraken_fil_df,fresa_kraken_fil_df$Treatment)

Healthy <-Treatment$healthy[c(1,2,3)]
Wilted <-Treatment$wilted[c(1,2,3)]

Healthy <- reshape(Healthy, idvar = "Sample",v.names= c("Abundance"), timevar = "OTU", direction = "wide")
rownames(Healthy) <- Healthy$Sample
Healthy <- select(Healthy, -Sample)
Wilted <- reshape(Wilted, idvar = "Sample",v.names= c("Abundance"), timevar = "OTU", direction = "wide")
rownames(Wilted) <- Wilted$Sample
Wilted <- select(Wilted, -Sample)

colH<-specaccum(Healthy, method = "collector")
plot(colH, xlab="Muestras", ylab="Abundancia", col="blue",ylim = c(8500,9000))
plotH + points(colH$richness, pch=19, col="darkblue")

colW<-specaccum(Wilted, method = "collector")
plot(colW, xlab="Muestras", ylab="Abundancia", col="blue",ylim = c(8500,9000))
points(colW$richness, pch=19, col="darkblue")

## Rarefaccion basada en individuos
NH <- colSums(Healthy)
subs3H <- c(seq(500, 4000, by = 500), sum(NH)) 
rar3H <- rarefy(NH, sample = subs3H, se = T, MARG = 2)
rar3H

NW <- colSums(Wilted)
subs3W <- c(seq(500, 4000, by = 500), sum(NW)) 
rar3W <- rarefy(NW, sample = subs3W, se = T, MARG = 2)
rar3W










################################################################################

# https://blogs.iu.edu/ncgas/2019/09/04/rarefaction-curves-for-metagenomic-datasets/ 

#Ejemplo rarefaccion con BCI 
#https://ciespinosa.github.io/AlphaDiversidad/rarefaccion.html
data(BCI)
dim(BCI)
# subconjunto de BCI tomando 10 muestras aleatorias
BCI_sub <- BCI[c(sample(1:50, 10, replace = TRUE)),]
dim(BCI_sub)
# eliminamos columnas de otus con suma de conteos en ceros
BCI_sub <- BCI_sub[,colSums(BCI_sub)>=1 ] 
dim(BCI_sub)
# realizamos una curva de acumulación de especies
col<-specaccum(BCI_sub, method = "collector")
plot(col, xlab="Parcelas", ylab="Número de especies", col="blue")
points(col$richness, pch=19, col="darkblue")
# "aleatorizar el orden y construir una curva, este proceso es conocido como rarefacción"
# La función rarefy arroja como resultado la riqueza de especies esperada en un determinado tamaño de muestra.
# la rarefacción permite hacer una interpolación de los datos, obteniendo una riqueza esperada en un tamaño de 
# muestra menor al tamaño que hemos logrado, de esta forma este proceso nos da no solamente la riqueza sino un error estándar.

## Rarefaccion basada en individuos
# Sumamos la abundancia de cada especie
N <- colSums(BCI_sub)
# Hacemos un vector con los tamaños de muestra sobre los cuales haremos 
# la interpolación. El dato final de este vector es el tamaño total de la 
# muestra (sum(N))
subs3 <- c(seq(500, 4000, by = 500), sum(N)) 
# Ejecutamos la rarefacción
rar3 <- rarefy(N, sample = subs3, se = T, MARG = 2)
rar3
# muestra la cantidad deespecies que se espera tener para diferentes tamaños de muestras

## Rarefaccion basada en muestras
#acomulacion de especies
rand<-specaccum(BCI_sub, method = "random", permutations=100)
rand
#vemos la acomulacion de otus encontrados dependiendo la cantidad de muestras
