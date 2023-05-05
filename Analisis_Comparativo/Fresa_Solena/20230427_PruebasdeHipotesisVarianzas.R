## PRUEBA DE HIPOTESIS - ANALISIS DE VARIANZAS PARA DIVERSIDADES CHAO1

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

## Calculamos la diversidad Chao1
Chao1_OTU <- estimateR(t(OTU))  

Chao1_OTU_df <- data.frame(sample=colnames(Chao1_OTU),value=Chao1_OTU[2, ])#, measure=rep('Chao1',length(Chao1_OTU)))

total <-cbind(Chao1_OTU_df,SAM)

# media por grupos
mu <- ddply(total, "Treatment", summarise, grp.mean=mean(value))
# numero de muestras
n <- total%>%count('Treatment')
# varianzas 
sigma <- ddply(total, "Treatment", summarise, s.var=var(value))

# Estadistico de prueba (F de Fisher)
F <- (sigma[1,2])/(sigma[2,2])

# grados de liberad
gl <- (n[1,2]-1)/(n[2,2]-1)

p<-ggplot(total, aes(x=value))+
  geom_histogram(color="pink",fill="black")+
  facet_grid(Treatment ~ .)
p+geom_vline(data=sigma, aes(xintercept=s.var, color="red"),
             linetype="dashed")

ggplot(total, aes(x=Treatment, y=value, color=Treatment)) + geom_point(size=2)






#################################################################################
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

## ACTINOBACTERIA
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


##AQUI QUEDE.......20230404_20:42---------------------------------------------------------------------
## Calculamos la diversidad Chao1
Chao1_OTU <- estimateR(glom_df2$Abundance)  

Chao1_OTU_df <- data.frame(sample=colnames(Chao1_OTU),value=Chao1_OTU[2, ])#, measure=rep('Chao1',length(Chao1_OTU)))

total <-cbind(Chao1_OTU_df,SAM)

# media por grupos
mu <- ddply(total, "Treatment", summarise, grp.mean=mean(value))
# numero de muestras
n <- total%>%count('Treatment')
# varianzas 
sigma <- ddply(total, "Treatment", summarise, s.var=var(value))

# Estadistico de prueba (F de Fisher)
F <- (sigma[1,2])/(sigma[2,2])

# grados de liberad
gl <- (n[1,2]-1)/(n[2,2]-1)

p<-ggplot(total, aes(x=value))+
  geom_histogram(color="pink",fill="black")+
  facet_grid(Treatment ~ .)
p+geom_vline(data=sigma, aes(xintercept=s.var, color="red"),
             linetype="dashed")


