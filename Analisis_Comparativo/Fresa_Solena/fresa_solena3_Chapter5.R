## ANALISIS DE POTENCIA - PRUEBA DE HIPOTESIS

library("phyloseq")
library("ggplot2")
library("vegan")
library("RColorBrewer")
library("dplyr")
library("plyr")

##cargamos la tabla de abundancia
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

OTU <- fresa_kraken_fil@otu_table@.Data
SAM <- fresa_kraken_fil@sam_data
OTU_SAM <-merge(OTU,SAM)
OTU_SAM <-t(OTU_SAM)
## UNIR TABLA DE ABUNDANCIAS CON METADATA DE SANOS Y ENFERMOS

## Calculamos la diversidad Shannon 

## Se usa la funcion diversidad del paquete vegan para calcular el indice Shannon
## Se realiza el dataframe del indice de Shannon
Shannon_OTU_SAM <- diversity(OTU_SAM, "shannon")
Shannon_OTU_df <- data.frame(sample=names(Shannon_OTU),value=Shannon_OTU,measure=rep("Shannon", length(Shannon_OTU)))





#dividir la trama en múltiples paneles
p<-ggplot(df_H_G93BUm3, aes(x=value))+
  geom_histogram(color="black",fill="black")+
  facet_grid(Group ~ .)
ggsave("paneles.png", plot = p)

#Calcular la media de cada grupo
#calcular la diversidad promedio de Shannon de cada grupo usando el paquete Plyr
##  GRUPOS SANOS VS. ENFERMOS 
mu <- ddply(df_H_G93BUm3, "Group", summarise, grp.mean=mean(value))
head(mu)
#agrega las lineas de la media
p2<-p+geom_vline(data=mu, aes(xintercept=grp.mean, color="red"),
             linetype="dashed")

ggsave("paneles.png", plot = p2)

##5.2.3 Cálculo de la potencia o el tamaño de la muestra mediante la función R power.t.test
#Dado que se desconoce la desviación estándar de la diferencia de medias,
#es necesario estimarlo

aggregate(value ~ Group, 
          df_H_G93BUm3,
          mean)

#Group value
#1   BUm3to3.5 2.504
#2 NOBUm3to3.5 2.205

aggregate(value ~ Group, df_H_G93BUm3,var)

#Group   value
#1   BUm3to3.5 0.02892
#2 NOBUm3to3.5 0.04349

n1 <- 9
n2 <-7
s1<-sqrt(0.02892) # raiz cuadrada de lavarianza 
s1 # 'desviacion estandar'
#[1] 0.1701
s2<-sqrt(0.04349)
s2
#[1] 0.2085
s=sqrt((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2)
s # desviacion estandar tomando en cuenta ambas 'n1' y 'n2'
#[1] 0.05012

#prueba T-test
power.t.test(n=2:10,delta=2.504-2.205,sd=0.05012)
df_P <-data.frame(n,power)
df_P

#n = c(2, 3, 4, 5, 6, 7, 8, 9, 10)  
#power = c(0.8324, 0.9994, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000)

power <- sapply(n, function (x) power.t.test(n=x, delta=2.504-2.205,sd=0.05012)$power)
po<-plot(n, power, xlab  = "Sample Size per group", ylab  = "Power to reject null",
     main="Power curve for\n t-test with delta = 0.05",
     lwd=2, col="red", type="l")

abline(h = 0.90, col="blue")

ggsave("power.png",plot=po)

power.t.test(n=2:10,delta=2.504-2.205,sd=0.05012, type = "one.sample" )


#5.3 Power Analysis for Comparing Diversity Across More than Two Groups Using 
#    ANOVA	
df_H_G93WTm1N4 <- filter(df_H_G6,Group%in%c("G93m1","WTm1","G93m4","WTm4"))
df_H_G93WTm1N4 

fit = lm(formula = value~Group,data=df_H_G93WTm1N4)
anova (fit)


library(pwr)
pwr.anova.test(f= 0.23,k=4,n=45:55,sig.level=0.05)