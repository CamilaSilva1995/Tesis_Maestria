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
write.csv(total, "Medidas_Chao1.csv")

# media por grupos
mu <- ddply(total, "Treatment", summarise, grp.mean=mean(value))
# numero de muestras
n <- total%>%count('Treatment')
# varianzas 
sigma <- ddply(total, "Treatment", summarise, s.var=var(value))

# Estadistico de prueba (F de Fisher)
F <- (sigma[1,2])^2/(sigma[2,2])^2

# grados de liberad
gl <- (n[1,2]-1)/(n[2,2]-1)

p<-ggplot(total, aes(x=value))+
  geom_histogram(color="pink",fill="black")+
  facet_grid(Treatment ~ .)
p+geom_vline(data=sigma, aes(xintercept=s.var, color="red"),
             linetype="dashed")

ggplot(total, aes(x=Treatment, y=value, color=Treatment)) + geom_point(size=2)

alfa <- 0.05

# prueba de Fisher para igualdad de varianzas 
# la varianza mayor debe ser S1
s2<-sigma[1,2]
s1<-sigma[2,2]
n2<-n[1,2]
n1<-n[2,2]
# grados de liberad
v1<-n1-1
v2<-n2-1
gl <- v1+v2

#estadístico F calculado
Fs<-(s1)^2/(s2)^2
#F teórico para alfa = 0.05 (dos colas)
# c(0.025,0.975) secuencia de probabilidades 
Ftablal <- qf(c(alfa/2,alfa,2*alfa), v1, v2, lower.tail = TRUE)
Ftablar <- qf(c(alfa/2,alfa,2*alfa), v1, v2, lower.tail = FALSE)
#qf() function gives the quantile function
#function in R Language is used to compute the value of quantile function over F distribution for a sequence of numeric values. It also creates a density plot of quantile function over F Distribution.
Pl <- pf(Fs, v1, v2,lower.tail = TRUE)#[0,Fs]
Pr <- pf(Fs, v1, v2,lower.tail = FALSE)#[Fs,+inf]
#pf() function gives the distribution function
#We use the pf() to calculate the area under the curve for the interval [0,Fs] and [Fs,+inf]
Pl+Pr==1
#plot a histogram and compare it to the probability density function of the F-distribution with v1 and v2 (pink line).
x <- rf(100000, v1, v2)
hist(x, 
     breaks = 'Scott', 
     freq = FALSE, 
     xlim = c(0,3), 
     ylim = c(0,1),
     xlab = '')
curve(df(x, v1, v2), from = 0, to = 4, n = 5000, col= 'pink', lwd=2, add = T)

##EJEMPLO
n1 <- 10
n2 <- 20
s1_2<-26.4
s2_2<-12.7
v1<-n1-1
v2<-n2-1
gl <- v1+v2
alfa<-0.05
F <- 3#s1_2/s2_2
# Valores critico
Ftablal <- qf(c(alfa/2,alfa,2*alfa), v1, v2, lower.tail = TRUE)
Ftablar <- qf(c(alfa/2,alfa,2*alfa), v1, v2, lower.tail = FALSE)
# Regiones de rechazo 
Pl <- pf(F, v1, v2,lower.tail = TRUE)#[0,Fs]
Pr <- pf(F, v1, v2,lower.tail = FALSE)#[Fs,+inf]
Pl+Pr==1
x <- rf(100000, v1, v2)
hist(x, 
     breaks = 'Scott', 
     freq = FALSE, 
     xlim = c(0,3), 
     ylim = c(0,1),
     xlab = '')
curve(df(x, v1, v2), from = 0, to = 4, n = 5000, col= 'pink', lwd=2, add = T)

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
Data <- glomToGraph(merge_Bacteria,'Genus')
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


##AQUI QUEDE.......
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


