## PRUEBA DE HIPOTESIS SOBRE MEDIAS (GENERAL Y PARA LOS TAXONES ENFOQUE)
#https://fhernanb.github.io/Manual-de-R/ph.html

library("phyloseq")
library("ggplot2")
library("vegan")
library("RColorBrewer")
library("stringi")
library("dplyr")
library("plyr")
library("gginference")

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

## Prueba de hipotesis con medias con indice Shannon

fresa_kraken_df <- psmelt(fresa_kraken_fil)

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
write.csv(total, "Medidas_Shannon_total.csv")

# media por grupos
mu <- ddply(total, "Treatment", summarise, grp.mean=mean(value))

p<-ggplot(total, aes(x=value))+
  geom_histogram(color="pink",fill="black")+
  facet_grid(Treatment ~ .)
p+geom_vline(data=mu, aes(xintercept=grp.mean, color="red"),
             linetype="dashed")
ggsave("medias_Shannon.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

# se fija el nivel de significancia
alfa <- 0.05
# numero de muestras
n <- total%>%count('Treatment')
# varianzas 
sigma <- ddply(total, "Treatment", summarise, s.var=var(value))

# S_{p} = \sqrt{\frac{(n_{1}-1)S_{1}^{2} + (n_{2}-1)S_{2}^{2}}{n_{1} + n_{2} - 2}}
#Sp <-  sqrt(((n[1,2]-1)*(sigma[1,2]*sigma[1,2]) + (n[2,2]-1)*(sigma[2,2]*sigma[2,2]))/(n[1,2]+n[2,2]-2)) 
n1<-n[1,2]
n2<-n[2,2]
# grados de liberad
v1<-n1-1
v2<-n2-1
gl <- v1+v2
# varianza muestral (sigma)
ss1<-sigma[1,2]*v1
ss2<-sigma[2,2]*v2

Sp2 <- (ss1+ss2)/(v1+v2)

# Estadistico de prueba (t-Student)
# T = \frac{\bar{Y_{1}} - \bar{Y_{2}}}{S_{p}\sqrt{\frac{1}{n_{1}} + \frac{1}{n_{2}}}}
T <- (mu[1,2] - mu[2,2])/ (sqrt(Sp2/n1 + Sp2/n2))

# P-value -> probabilidad para rechazar la hipotesis (nivel de significancia observada)

p_value<-2*pt(q=T,df=v1+v2,lower.tail = FALSE)


### no se rechaza H0 por que p-value>alfa
# esto es suponiendo varianzas iguales

totalH <- total[total$Treatment == "healthy", ]
totalW <- total[total$Treatment == "wilted", ]

pruebat <- t.test(totalH$value, totalW$value, var.equal = TRUE, alternative = "two.sided")
pruebat
ggttest(pruebat)
ggsave("pruebat_varianzasiguales_Shannon.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

### suponiendo varianzas diferentes

pruebat2 <- t.test(totalH$value, totalW$value, var.equal = FALSE, alternative = "two.sided")
pruebat2
ggttest(pruebat2)
ggsave("pruebat_varianzasdiferentes_Shannon.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

####

qqnorm(totalH$value,main = "Healthy");qqline(totalH$value)
qqnorm(totalW$value,main = "Wilted");qqline(totalW$value)

#con shapiro para ver la distribucion de nuestros datos
shapiro.test(totalH$value)
shapiro.test(totalW$value)

wilcox.test(totalH$value, totalW$value)

datos <- sort(total$value)
rango <- rank(datos)
observaciones <- cbind(total$value,datos,rango)

##############################
# prueba de Fisher para igualdad de varianzas 
s1<-sigma[1,2]
s2<-sigma[2,2]

F=s1/s2

Ftabla <- qf(c(0.025,0.975), v1, v2)

Pr <- pf(F, v1, v2)

################################
###https://www.jmp.com/es_co/statistics-knowledge-portal/t-test/two-sample-t-test.html
###https://rstudio-pubs-static.s3.amazonaws.com/423172_fbe936199b504e38882f26ececc66a2a.html#prueba-t-con-varianzas-diferentes
nH <- totalH%>%count('Treatment')
nH <- nH[1,2]
nW <- totalW%>%count('Treatment')
nW <- nW[1,2]
P1 = 1 - alfa
P2 = 1 - alfa/2

t1 <- qt(P1, gl)
t2 <- qt(P2, gl)

#The function qt returns the value of the inverse cumulative density function (cdf) of the Student t distribution given a certain random variable x and degrees of freedom df. 
# Una cola
sprintf("t_1cola: %g", t1)
# Dos colas
sprintf("t_2colas: %g", t2)

# p-errorI para la t-calculada
Ptcalc1 <- pt(T, gl, lower.tail = FALSE)

# The function pt returns the value of the cumulative density function (cdf) of the Student t distribution given a certain random variable x and degrees of freedom df.
sprintf("P_errorI_1cola: %g", Ptcalc1)
sprintf("P_errorI_2cola: %g", Ptcalc1*2)

ggplot(total, aes(x=Treatment, y=value, color=Treatment)) + geom_point(size=2)

#################################################################################
## TAXONES ENFOQUE
## FUSARIUM

## Subconjunto de "Eukaryota"
merge_Eukaryota<-subset_taxa(fresa_kraken_fil,Kingdom=="Eukaryota")

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
Shannon_OTU_df <- data.frame(sample=names(Shannon_OTU),value=Shannon_OTU,measure=rep("Shannon", length(Shannon_OTU)))

Simp_OTU <- diversity(df, "simpson")
Simp_OTU_df <- data.frame(Simpson=Simp_OTU)

#unir con 
total_Shannon <-cbind(glom@sam_data,Shannon_OTU_df)
total <-cbind(glom@sam_data,Shannon_OTU_df,Simp_OTU_df )
write.csv(total_Shannon, "Medidas_Shannon_Genero.csv")

# media por grupos para indice Shannon
mu_Shannon <- ddply(total_Shannon, "Treatment", summarise, grp.mean=mean(value))
# numero de muestras
n <- total_Shannon%>%count('Treatment')
# varianzas 
sigma <- ddply(total_Shannon, "Treatment", summarise, s.var=var(value))

# S_{p} = \sqrt{\frac{(n_{1}-1)S_{1}^{2} + (n_{2}-1)S_{2}^{2}}{n_{1} + n_{2} - 2}}
#Sp <-  sqrt(((n[1,2]-1)*(sigma[1,2]*sigma[1,2]) + (n[2,2]-1)*(sigma[2,2]*sigma[2,2]))/(n[1,2]+n[2,2]-2)) 
n1<-n[1,2]
n2<-n[2,2]
# grados de liberad
v1<-n1-1
v2<-n2-1
gl <- v1+v2
# varianza muestral (sigma)
ss1<-sigma[1,2]*v1
ss2<-sigma[2,2]*v2

Sp2 <- (ss1+ss2)/(v1+v2)

# Estadistico de prueba (t-Student)
# T = \frac{\bar{Y_{1}} - \bar{Y_{2}}}{S_{p}\sqrt{\frac{1}{n_{1}} + \frac{1}{n_{2}}}}
T <- (mu_Shannon[1,2] - mu_Shannon[2,2])/ (sqrt(Sp2/n1 + Sp2/n2))


# # S_{p} = \sqrt{\frac{(n_{1}-1)S_{1}^{2} + (n_{2}-1)S_{2}^{2}}{n_{1} + n_{2} - 2}}
# Sp <-  sqrt(((n[1,2]-1)*(sigma[1,2]*sigma[1,2]) + (n[2,2]-1)*(sigma[2,2]*sigma[2,2]))/(n[1,2]+n[2,2]-2)) 
# 
# # Estadistico de prueba (t-Student)
# # T = \frac{\bar{Y_{1}} - \bar{Y_{2}}}{S_{p}\sqrt{\frac{1}{n_{1}} + \frac{1}{n_{2}}}}
# T <- (mu_Shannon[1,2] - mu_Shannon[2,2])/ (Sp * sqrt(1/n[1,2] + 1/n[2,2]))

# nivel de significancia (alfa=0.05)
# P-value -> probabilidad para rechazar la hipotesis

p <- ggplot(total_Shannon, aes(x=value))+
  geom_histogram(color="pink",fill="black")+
  facet_grid(Treatment ~ .)
q <- p + geom_vline(data=mu_Shannon, aes(xintercept=grp.mean, color="red"),linetype="dashed")
ggsave("mediasFusarium_Shannon.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")


total_ShannonH <- total_Shannon[total_Shannon$Treatment == "healthy", ]
total_ShannonW <- total_Shannon[total_Shannon$Treatment == "wilted", ]


# esto essuponiendo varianzas iguales

pruebat <- t.test(total_ShannonH$value, total_ShannonW$value, var.equal = TRUE, alternative = "two.sided")
pruebat
ggttest(pruebat)
ggsave("pruebat_varianzasdiferentesFusarium_Shannon.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")


### suponiendo varianzas diferentes

pruebat2 <- t.test(total_ShannonH$value, total_ShannonW$value, var.equal = FALSE, alternative = "two.sided")
pruebat2
ggttest(pruebat2)
ggsave("pruebat_varianzasdiferentesFusarium_Shannon.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

###

P1 = 1 - alfa
P2 = 1 - alfa/2

t1 <- qt(P1, gl)
t2 <- qt(P2, gl)
#The function qt returns the value of the inverse cumulative density function (cdf) of the Student t distribution given a certain random variable x and degrees of freedom df. 

# Una cola
sprintf("t_1cola: %g", t1)
# Dos colas
sprintf("t_2colas: %g", t2)

# p-errorI para la t-calculada
Ptcalc1 <- pt(T, gl, lower.tail = FALSE)
# The function pt returns the value of the cumulative density function (cdf) of the Student t distribution given a certain random variable x and degrees of freedom df.
sprintf("P_errorI_1cola: %g", Ptcalc1)
sprintf("P_errorI_2cola: %g", Ptcalc1*2)

library('gginference')
ggttest(pruebat)





mu_Simp <- ddply(total, "Treatment", summarise, grp.mean=mean(Simpson))

p2 <- ggplot(total, aes(x=Simpson))+
  geom_histogram(color="pink",fill="black")+
  facet_grid(Treatment ~ .)
q2 <- p2 + geom_vline(data=mu_Simp, aes(xintercept=grp.mean, color="red"),linetype="dashed")


ggplot(total, aes(x=Treatment, y=Shannon, color=Treatment)) + geom_point(size=2)
