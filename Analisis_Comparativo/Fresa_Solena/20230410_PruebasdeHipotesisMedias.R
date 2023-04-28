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


OTU <- fresa_kraken_fil@otu_table@.Data
SAM <- fresa_kraken_fil@sam_data

## UNIR TABLA DE ABUNDANCIAS CON METADATA DE SANOS Y ENFERMOS

## Calculamos la diversidad Shannon 

## Se usa la funcion diversidad del paquete vegan para calcular el indice Shannon
## Se realiza el dataframe del indice de Shannon
OTU <- t(OTU)
Shannon_OTU <- diversity(OTU, "shannon")
Shannon_OTU_df <- data.frame(sample=names(Shannon_OTU),value=Shannon_OTU,measure=rep("Shannon", length(Shannon_OTU)))
total <-cbind(Shannon_OTU_df,SAM)

# media por grupos
mu <- ddply(total, "Treatment", summarise, grp.mean=mean(value))

p<-ggplot(total, aes(x=value))+
  geom_histogram(color="pink",fill="black")+
  facet_grid(Treatment ~ .)
p+geom_vline(data=mu, aes(xintercept=grp.mean, color="red"),
             linetype="dashed")

ggplot(total, aes(x=Treatment, y=value, color=Treatment)) + geom_point(size=2)


###EJEMPLO 10.14..PAG 524 DEL LIBRO DE ESTADISTICA
# Prueba de medias para Shannon



#################################################################################

## Subconjunto de "Eukaryota"
merge_Eukaryota<-subset_taxa(fresa_kraken_fil,Kingdom=="Eukaryota")
## Subconjunto de "Bacteria"
merge_Bacteria<-subset_taxa(fresa_kraken_fil,Kingdom=="Bacteria")

glomToGraph<-function(phy,tax){
  ## creamos el subconjunto dependiendo del linaje taxonomico deseado
  glom <- tax_glom(phy, taxrank = tax)
  glom_df <- psmelt(glom)
  ## sacamos los porcentajes
  percentages <- transform_sample_counts(glom, function(x) x*100 / sum(x) )
  percentages_df <- psmelt(percentages)
  return(list(glom,glom_df,percentages,percentages_df))
}

# PRUEBAS A NIVEL DE GENERO
Data <- glomToGraph(merge_Eukaryota,'Genus')
glom <- Data[[1]] # phyloseq
glom_df <- Data[[2]] # dataframe
percentages <- Data[[3]] # phyloseq
percentages_df <- Data[[4]] # dataframe

## solo tomamos las columnas Sample, OTU y Abundance
glom_df2 <- glom_df[c(1,2,3)]
head(glom_df2)

# queremos pasar de dataframe a table, para calcular el alfa diversidad
df <- reshape(glom_df2, idvar = "Sample",v.names= c("Abundance"), timevar = "OTU", direction = "wide")
rownames(df) <- df$Sample
df <- select(df, -Sample)
head(df)

## Calculamos la diversidad Shannon 
## Se usa la funcion diversidad del paquete vegan para calcular el indice Shannon
## Se realiza el dataframe del indice de Shannon

Chao1_OTU <- estimateR(glom_df2$Abundance)  

Shannon_OTU <- diversity(df, "shannon")
Shannon_OTU_df <- data.frame(Shannon=Shannon_OTU)

Simp_OTU <- diversity(df, "simpson")
Simp_OTU_df <- data.frame(Simpson=Simp_OTU)

#unir con 
total <-cbind(glom@sam_data,Shannon_OTU_df,Simp_OTU_df )

mu_Shannon <- ddply(total, "Treatment", summarise, grp.mean=mean(Shannon))

p <- ggplot(total, aes(x=Shannon))+
  geom_histogram(color="pink",fill="black")+
  facet_grid(Treatment ~ .)
q <- p + geom_vline(data=mu_Shannon, aes(xintercept=grp.mean, color="red"),linetype="dashed")

mu_Simp <- ddply(total, "Treatment", summarise, grp.mean=mean(Simpson))

p2 <- ggplot(total, aes(x=Simpson))+
  geom_histogram(color="pink",fill="black")+
  facet_grid(Treatment ~ .)
q2 <- p2 + geom_vline(data=mu_Simp, aes(xintercept=grp.mean, color="red"),linetype="dashed")


ggplot(total, aes(x=Treatment, y=Shannon, color=Treatment)) + geom_point(size=2)


