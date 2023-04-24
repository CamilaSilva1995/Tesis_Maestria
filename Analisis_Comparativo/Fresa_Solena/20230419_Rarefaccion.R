library(vegan)


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
## "aleatorizar el orden y construir una curva, este proceso es conocido como rarefacción"
# La función rarefy arroja como resultado la riqueza de especies esperada en un determinado tamaño de muestra.
# la rarefacción permite hacer una interpolación de los datos, obteniendo una riqueza esperada en un tamaño de 
# muestra menor al tamaño que hemos logrado, de esta forma este proceso nos da no solamente la riqueza sino un error estándar.

## Rarefaccion basada en individuos
#Sumamos la abundancia de cada especie
N <- colSums(BCI_sub)
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
rand<-specaccum(BCI_sub, method = "random", permutations=100)
rand
#vemos la acomulacion de otus encontrados dependiendo la cantidad de muestras

#####
#rarefacción basada en individuos
#utilizando el valor de abundancia de la parcela con menos individuos. 
R.rar <- rarefy(BCI_sub, min(rowSums(BCI_sub)), se=TRUE)
R.rar #riqueza interpolada (con rarefacción) de cada una de las parcelas
#gráfico en el cual podamos comparar la riqueza total de las parcelas del BCI
par(mfcol=c(2,1), oma=c(2,1,1,1), mar=c(3,3,1,1))
Rt<- specnumber(BCI_sub)
barplot(Rt, ylim=c(0,120), col="grey30", cex.main=0.8,
        cex.names=0.7, cex.axis=0.7, mgp=c(2, 0.35, 0))
mtext("a) Riqueza total por parcela", side=3, cex=0.8)
mtext("Número de especies", side =2, at=0, line=2, cex=0.8)
x<-(barplot(R.rar[1,], ylim=c(0,120), col="grey30", 
            cex.main=0.8, cex.axis=0.7, cex.names=0.7, mgp=c(2, 0.35, 0), 
            xlab="Sitios", cex.lab=0.8))
mtext("b) Riqueza con rarefacción basada en individuos por parcela", side =3, cex=0.8)

#al comparar la riqueza entre las diferentes parcelas el patrón se mantiene entre los datos reales y los datos de la rarefacción



####
#La curva de rarefacción da la riqueza de especies anotadas es un gráfico, el
#número total de anotaciones de especies distintas en función del número de
#secuencias muestreadas.


library("phyloseq")
library("ggplot2")
library("vegan")
library("RColorBrewer")
library("stringi")
library("dplyr")
library("plyr")
library("edgeR")

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
outpath = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"
### Cargado de datos 
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










