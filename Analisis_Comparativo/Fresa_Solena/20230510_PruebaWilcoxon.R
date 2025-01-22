## PRUEBA DE WIlCOXON

library("phyloseq")
library("ggplot2")
library("vegan")
library("RColorBrewer")
library("stringi")
library("dplyr")
library("plyr")

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
outpath = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"
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
percentages_fil <- transform_sample_counts(fresa_kraken_fil, function(x) x*100 / sum(x) )
percentages_df <- psmelt(percentages_fil)

###EJEMPLO 10.16..PAG 532 DEL LIBRO DE ESTADISTICA
# Prueba de varianzas para Chao1

OTU <- fresa_kraken_fil@otu_table@.Data
SAM <- fresa_kraken_fil@sam_data

## UNIR TABLA DE ABUNDANCIAS CON METADATA DE SANOS Y ENFERMOS

#################################################
## UNIR TABLA DE ABUNDANCIAS CON METADATA DE SANOS Y ENFERMOS

## Calculamos la diversidad Shannon 

## Se usa la funcion diversidad del paquete vegan para calcular el indice Shannon
## Se realiza el dataframe del indice de Shannon

## WILCOX PARA SHANNON

OTU <- t(OTU)
Shannon_OTU <- diversity(OTU, "shannon")
Shannon_OTU_df <- data.frame(sample=names(Shannon_OTU),value=Shannon_OTU,measure=rep("Shannon", length(Shannon_OTU)))
total <-cbind(Shannon_OTU_df,SAM)

total$Rank<-rank(total$value)

ggplot(data = total, aes(x = Rank, y = value)) +
  geom_point(aes(colour = Treatment), size = 3) +
  ylab("") + xlab("rango") +
  theme_bw() +
  theme(axis.text.y = element_blank()) + 
  ggtitle("Muestras procedentes de la misma población")
ggsave("Wilcoxon_Shannon.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

healthy<-total[total$Treatment=="healthy",]$value
wilted<-total[total$Treatment=="wilted",]$value


wilcox.test(x = healthy, y = wilted, alternative = "two.sided", mu = 0,
            paired = FALSE, conf.int = 0.95)


## WILCOX PARA CHAO1

