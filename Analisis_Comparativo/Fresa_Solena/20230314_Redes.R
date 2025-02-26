### REDES
### Redes simples
### https://haydeeperuyero.github.io/Analisis_datos_microbioma_Eq4/Chapter7/Chapter7-html
### Redes de Coocurrencia
### http://www.castrolab.org/isme/microbial_networks/microbial_networks.html#workshop
### http://www.castrolab.org/isme/microbial_networks/microbial_networks.html#comau
### http://www.aic.uva.es/cuentapalabras/colocaci%C3%B3n-coocurrencia-y-redes-l%C3%A9xicas.html


library("phyloseq")
library("ggplot2")
library("igraph")
library("vegan")
library("RColorBrewer")
library("GUniFrac")

setwd("/home/camila/GIT/Tesis_Maestria/Data2/fresa_solena")
fresa_kraken <- import_biom("fresa_kraken.biom")
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data,4,100)
colnames(fresa_kraken@otu_table@.Data) <- substr(colnames(fresa_kraken@otu_table@.Data),1,6)
metadata_fresa <- read.csv2("/home/camila/GIT/Tesis_Maestria/Data2/fresa_solena/metadata.csv",header =  FALSE, row.names = 1, sep = ",")
fresa_kraken@sam_data <- sample_data(metadata_fresa)
fresa_kraken@sam_data$Sample<-row.names(fresa_kraken@sam_data)
colnames(fresa_kraken@sam_data)<-c('Treatment','Samples')
samples_to_remove <- c("MP2079","MP2080","MP2088","MP2109","MP2137")
fresa_kraken_fil <- prune_samples(!(sample_names(fresa_kraken) %in% samples_to_remove), fresa_kraken)
percentages_fil <- transform_sample_counts(fresa_kraken_fil, function(x) x*100 / sum(x) )
## podemos hacer un data.frame uniendo toda la informacion del objeto phyloseq
df <- psmelt(fresa_kraken_fil)

## Hay dos funciones en el paquete phyloseq para trazar la red del microbioma usando “ggplot2”: plot_network() y plot_net()

## Se crean un grafo basado en igraph, basado en el método de distancia por defecto, Jaccard y una distancia máxima entre nodos conectados de 0,8. El “Treatment” se utiliza para los mapeos de color y forma para visualizar la estructura de las muestras

## Hacemos una gráfica a partir de el objeto phyloseq
ig <- make_network(fresa_kraken_fil, max.dist=0.8)
## Graficamos
plot_network(ig, fresa_kraken_fil, color="Treatment", shape="Treatment")
 
## Los siguientes códigos crean una red basada en una distancia máxima entre los nodos conectados de 0,5
## Graficamos sin necesidad de hacer el objeto anterior
plot_net(fresa_kraken_fil, maxdist = 0.5, color = "Treatment", shape="Treatment")


################################################################################
## REDES DE COOCURRENCIA

#setwd("~/GIT/Tesis_Maestria/Data/fresa_solena/Data1/Redes")
#20230314-1_raw_network -> (correlacion normalizada de -1 a 1, con fusarion como taxon clave)

# crear exel de correlacion para 5506 - Fusarium, con los siguientes taxa:

# 5125 - Hypocreales
# 57161 - Fusarium decemcellulare
# 37994 - Beauveria felina
# 5141 - Neurospora crassa

#ejemplo de los datos de correlacion para cada genero correlacionado con fusarium
Fusarium <- c(1,2,3)

Hypocreales <- c(4,5,6)

cor(Fusarium,Hypocreales)

plot(Fusarium,Hypocreales)

# se puedem hacer pruebas de permutación y luego volver a calcular la correlacion
Hypocreales_p <- c(6,5,4)

cor(Fusarium,Hypocreales_p)

# al desordenar los puntos se vuelve a calcular la correlacion, y si la muestra no es importante se pueden permutar las columnas y no ver cambios notorios

# Prueba de Boostrap
# Desordena aleatoreamente y se puede repetir el calculo de la correlacion







### Gephi
# se pueden vizualizar los resultados con gephi

