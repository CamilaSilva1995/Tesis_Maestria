## ANALISIS DE POTENCIA - PRUEBA DE HIPOTESIS

### https://fhernanb.github.io/Manual-de-R/ph.html

library("phyloseq")
library("ggplot2")
library("vegan")
library("RColorBrewer")
library("dplyr")
library("plyr")

##cargamos la tabla de abundancia
setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena")
fresa_kraken <- import_biom("fresa_kraken.biom")
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data,4,100)
colnames(fresa_kraken@otu_table@.Data) <- substr(colnames(fresa_kraken@otu_table@.Data),1,6)
metadata_fresa <- read.csv2("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/metadata.csv",header =  FALSE, row.names = 1, sep = ",")
fresa_kraken@sam_data <- sample_data(metadata_fresa)
fresa_kraken@sam_data$Sample<-row.names(fresa_kraken@sam_data)
colnames(fresa_kraken@sam_data)<-c('Treatment','Samples')
samples_to_remove <- c("MP2079","MP2080","MP2088","MP2109","MP2137")
fresa_kraken_fil <- prune_samples(!(sample_names(fresa_kraken) %in% samples_to_remove), fresa_kraken)
percentages_fil <- transform_sample_counts(fresa_kraken_fil, function(x) x*100 / sum(x) )


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

#mesiapor grupos
mu <- ddply(total, "Treatment", summarise, grp.mean=mean(value))

p<-ggplot(total, aes(x=value))+
  geom_histogram(color="pink",fill="black")+
  facet_grid(Treatment ~ .)
p+geom_vline(data=mu, aes(xintercept=grp.mean, color="red"),
             linetype="dashed")


################################# CAP5 - POWER ################################# 

# ## Código tomado del libro Statistical analysis of microbiome data with R 
# ## Autores Yinglin Xia et al 
# ## Ejercicios transcrito por MArio Jardón y Janeth de Anda
# 
# ## Cargamos la tabla de abundancia
# abund_table = read.csv("ALSG93AGenus.csv",row.names=1, check.names=FALSE)
# abund_table_t<-t(abund_table)
# 
# ##calculamos la diversidad Shannon 
# library(vegan)
# #Se usa la funcion diversidad del paquete vegan para calcular el indice Shannon
# #Se realiza el dataframe del indice de Shannnon
# H<-diversity(abund_table_t, "shannon")
# df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon", length(H)))
# #Se obtiene la informacion agrupada para los datps de la muestra
# df_H$Group <- with(df_H, ifelse(as.factor(sample)%in% c("A11-28F","A12-28F",
#                                                         "A13-28F","A14-28F","A15-28F","A16-28F"),c("G93m1"), 
#                                 ifelse(as.factor(sample)%in% c("A21-28F","A22-28F","A23-28F","A24-28F",
#                                                                "A25-28F","A26-28F"),c("WTm1"),
#                                        ifelse(as.factor(sample)%in% c("C11-28F","C12-28F","C13-28F"),c("G93m4"),
#                                               ifelse(as.factor(sample)%in% c("C21-28F","C22-28F","C23-28F"),
#                                                      c("WTm4"),ifelse(as.factor(sample)%in% c("B11-28F",
#                                                                                               "B12-28F","B13-28F","B14-28F","B15-28F","D11-28F",
#                                                                                               "D12-28F","D13-28F","D14-28F"),c("BUm3to3.5"),
#                                                                       c("NOBUm3to3.5")))))))
# 
# df_H
# #El conjunto de datos incluye datos de muestra de los meses 1, 3, 3,5 y 4. 
# #Nos interesa en las comparaciones del tratamiento frente al no tratamiento 
# #durante los meses 3 y 3,5. Así que los datos de 3 a 3,5 meses como se indica a
# #continuación:
# library(dplyr)
# df_H_G6 <- select(df_H, Group,value)
# df_H_G93BUm3 <- filter(df_H_G6,Group=="BUm3to3.5"|Group=="NOBUm3to3.5")
# df_H_G93BUm3
# library(ggplot2)
# #dividir la trama en múltiples paneles
# p<-ggplot(df_H_G93BUm3, aes(x=value))+
#   geom_histogram(color="black",fill="black")+
#   facet_grid(Group ~ .)
# #Calcular la media de cada grupo
# #calcular la diversidad promedio de Shannon de cada grupo usando el paquete Plyr
# library(plyr)
# mu <- ddply(df_H_G93BUm3, "Group", summarise, grp.mean=mean(value))
# head(mu)
# #agrega las lineas de la media
# p+geom_vline(data=mu, aes(xintercept=grp.mean, color="red"),
#              linetype="dashed")
# 
# ##5.2.3 Cálculo de la potencia o el tamaño de la muestra mediante la función R power.t.test
# #Dado que se desconoce la desviación estándar de la diferencia de medias,
# #es necesario estimarlo
# 
# aggregate(value ~ Group, df_H_G93BUm3, mean)
# 
# #Group value
# #1   BUm3to3.5 2.504
# #2 NOBUm3to3.5 2.205
# 
# aggregate(value ~ Group, df_H_G93BUm3,var)
# 
# #Group   value
# #1   BUm3to3.5 0.02892
# #2 NOBUm3to3.5 0.04349
# 
# n1 <- 9
# n2 <-7
# s1<-sqrt(0.02892)
# s1
# #[1] 0.1701
# s2<-sqrt(0.04349)
# s2
# #[1] 0.2085
# s=sqrt((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2)
# s
# #[1] 0.05012
# 
# power.t.test(n=2:10,delta=2.504-2.205,sd=0.05012)
# df_P <-data.frame(n,power)
# df_P
# 
# n = c(2, 3, 4, 5, 6, 7, 8, 9, 10)  
# power = c(0.8324, 0.9994, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000)
# 
# power <- sapply(n, function (x) power.t.test(n=x, delta=2.504-2.205,sd=0.05012)$power)
# plot(n, power, xlab  = "Sample Size per group", ylab  = "Power to reject null",
#      main="Power curve for\n t-test with delta = 0.05",
#      lwd=2, col="red", type="l")
# abline(h = 0.90, col="blue")
# power.t.test(n=2:10,delta=2.504-2.205,sd=0.05012, type = "one.sample" )
# 
# #5.3 Power Analysis for Comparing Diversity Across More than Two Groups Using 
# #    ANOVA	
# df_H_G93WTm1N4 <- filter(df_H_G6,Group%in%c("G93m1","WTm1","G93m4","WTm4"))
# df_H_G93WTm1N4 
# 
# fit = lm(formula = value~Group,data=df_H_G93WTm1N4)
# anova (fit)
# 
# library(pwr)
# pwr.anova.test(f= 0.23,k=4,n=45:55,sig.level=0.05)
# 
# 
# ##Análisis de potencia para comparar un taxón de interés entre grupos. Función power.prop.test()
# 
# #Se carga una tabla con la abundancia de Butyrivibrio ﬁbrisolvens por muestra. 
# 
# abund_table_Spe=read.csv("https://raw.githubusercontent.com/MJardon913/Analisis_de_microbioma_en_R/main/Data/ALSG93AButyrivibrioSpecies.csv", row.names=1, check.names=FALSE)
# abund_table_Spe<-t(abund_table_Spe)
# 
# #Las muestras se agrupan según si tuvieron el tratamiento de Butirato o son del grupo de control.
# 
# grouping<-data.frame(row.names=rownames(abund_table_Spe),t(as.data.frame(strsplit(rownames(abund_table_Spe),"-"))))
# grouping$Group <- with(grouping,ifelse(as.factor(X1)%in%c("B11","B12","B13","B14","B15","D11","D12","D13","D14"),c("Butyrate"), c("Control")))
# Butyrivibrio_G <-cbind(abund_table_Spe, grouping)
# Butyrivibrio_G
# 
# #Se agrega una columna que clasifica las muestras según si Butyrivibrio está presente.
# 
# Butyrivibrio_G$Present <- ifelse((Butyrivibrio_G$Butyrivibrio > 0),"Present","Absent")
# Butyrivibrio_G
# 
# #Se crea una tabla de contingencia según los grupos de muestras y la presencia de Butyrivibrio.
# 
# library(MASS)
# tbl = table(Butyrivibrio_G$Group, Butyrivibrio_G$Present)
# tbl
# 
# #Se calcula la potencia para muestras de tamaño de 10 a 20, valor alfa = 0.05. En este caso p1 es la proporción de muestras con butirato que presentan Butirivibrio y p2 es la misma proporción entre las muestras de control.
# 
# power.prop.test(n=10:20, p1=1, p2=.57, sig.level=0.05, power=NULL, alternative=c("one.sided"), strict = FALSE)
# 
# ##Análisis de potencia usando el test \chi²
# 
# #Este análisis que se llevará a cabo con la función pwr.chisq.test() desde la tabla de contingencia necesita del cálculo del tamaño del efecto. Para esto se usará la función cramersV().
# 
# #install.packages("lsr")
# library(lsr)
# 
# #Cálculo del tamaño de efecto
# cramersV(tbl)
# 
# library(pwr)
# pwr.chisq.test(w = 0.3833108, N = 45:60, df = 1, sig.level = 0.05, power = NULL)
# 
# #También desde las proporciones 1 y .57 podemos hacer un análisis de potencia desde el test de Fisher.
# 
# install.packages("Exact")
# library(Exact)
# power.exact.test(1.0, 0.57, 15, 15, method="Fisher")
# 
# #Comparando la frecuencia de todos los taxa con el modelo Dirichlet-multinomial
# 
# #El paquete HMP contiene la función que será usada para estimar un modelo Dirichlet-multinomial desde nuestros datos.
# 
# #install.packages("HMP",repo="http://cran.r-project.org", dep=TRUE)
# library(HMP)
# 
# #Cargamos los datos con frecuencias por género de cinco muestras con tratamiento de Butirato y de cinco muestras del grupo control.
# 
# Buty=read.csv("https://raw.githubusercontent.com/MJardon913/Analisis_de_microbioma_en_R/main/Data/ALSG93A3.5mButyrateGenus.csv",row.names=1,check.names=FALSE)
# NOButy=read.csv("https://raw.githubusercontent.com/MJardon913/Analisis_de_microbioma_en_R/main/Data/ALSG93A3.5mNoButyrateGenus.csv",row.names=1, check.names=FALSE)
# head(Buty)
# head(NOButy)
# 
# #Usamos las función DM.MoM para obtener los parámetros de ambos modelos Dirchlet-multinomial.
# 
# Buty_t <- t(Buty)
# NOButy_t <- t(NOButy)
# fit_Buty <- DM.MoM(Buty_t);fit_NOButy <- DM.MoM(NOButy_t)
# 
# ##Análisis por composición de taxa
# 
# #Con la función MC.Xdc.statistics() se llevará a cabo este análisis. Dicha función necesita especificar la cantidad de experimentos MonteCarlo = numMC, y el número de reads por muestra = group_Nrs
# 
# numMC <- 1000
# nrsGrp1 <- rep(1000, 10)
# nrsGrp2 <- rep(1000, 10)
# group_Nrs <- list(nrsGrp1, nrsGrp2)
# 
# #De este modo se calcula el error tipo 1.
# 
# alphap <- fit_Buty$gamma
# pval1 <- MC.Xdc.statistics(group_Nrs, numMC, alphap, "hnull")
# pval1
# 
# #De este modo se calcula la potencia. 
# 
# alphap <- rbind(fit_Buty$gamma, fit_NOButy$gamma)
# pval2 <- MC.Xdc.statistics(group_Nrs, numMC, alphap)
# pval2
# 
# ##Análisis RAD
# 
# #En este caso se filtrarán las muestras con los diez taxa más abundantes, declarando al resto como "other".
# 
# filter_Buty<- Data.filter(Buty_t, "sample", 1000, 10)
# head(filter_Buty)
# 
# filter_NOButy<- Data.filter(NOButy_t, "sample", 1000, 10)
# head(filter_NOButy)
# 
# #Se entrenan modelos Dirichlet-multinomial con los datos filtrados
# 
# fit_Buty <- DM.MoM(filter_Buty);fit_NOButy <- DM.MoM(filter_NOButy)
# 
# #Del siguiente modo podemos ver los parámetros de estos modelos.
# 
# fit_Buty$pi
# fit_NOButy$pi
# fit_Buty$theta
# fit_NOButy$theta
# 
# pi0 <- fit_Buty$pi
# group_theta <- c(0.007523, 0.01615)
# 
# #De este modo podemos calcular el error tipo 1.
# 
# pval1 <- MC.Xmc.statistics(group_Nrs, numMC, pi0, group.theta=group_theta, type="hnull")
# pval1
# 
# #Así podemos calcular la potencia. 
# 
# group_pi <- rbind(fit_Buty$pi, fit_NOButy$pi)
# pval2 <- MC.Xmc.statistics(group_Nrs, numMC, pi0, group_pi, group_theta)
# pval2
# 
# #Otra función que se puede usar para estos cálculos es la función MC.Xmcupo.statistics.
# #Para el error tipo 1:
# 
# pval1 <- MC.Xmcupo.statistics(group_Nrs, numMC, pi0, group.theta=group_theta, type="hnull")
# pval1
# 
# #Para la potencia:
# 
# pval2 <- MC.Xmcupo.statistics(group_Nrs, numMC, group.pi=group_pi, group.theta=group_theta)
# pval2
# 
# ##Cálculo del tamaño del efecto
# 
# group_data <- list(filter_Buty, filter_NOButy)
# effect <- Xmcupo.effectsize(group_data)
# effect
# 


# 
# #-----------------------------Medidas y cálculos para la diversidad alfa
# 
# #Primero, leamos los datos de abundacia del conteo de generos de R y cargaremos el paquete Vegan.
# options(width=65,digits=4)
# abund_table=read.csv("https://raw.githubusercontent.com/MJardon913/Analisis_de_microbioma_en_R/main/Data/ALSG93A3.5mButyrateGenus.csv",row.names=1,check.names=FALSE)
# library(vegan)
# 
# #El siguiente código de pocas líneas muestra que la estructura de datos is taxa (en este caso, géneros) por formato de muestra.
# 
# head(abund_table)
# 
# #La tabla de datos debe transformarse en muestras por taxones (géneros) antes de calcular varias diversidades.
# abund_table<-t(abund_table)
# head(abund_table)
# 
# #Aquí se usa la función specnumber() para calcular el número de géneros
# 
# num_genera<-specnumber(abund_table)
# num_genera
# 
# #--------------------------------Índice Chao1
# 
# #Ahora estimaremos el número de géneros y el índice Chao1 en nuestro conjunto de datos, usando la función estimate()
# index=estimateR(abund_table)
# index
# 
# #La función anterior genera 5 índices, del cual tomamos Chao1 usando el siguiente código en R
# chao1_genus=estimateR(abund_table)[2,]
# chao1_genus
# 
# #---------------------------------Índice de Shannon-Wiener
# 
# #Usamos la función de diversidad en el paquete Vegan para calcular el índice de Shannon-Wiener, para ello tenemos el siguiente código en R:
# shannon_genus<-diversity(abund_table, index = "shannon", MARGIN=1)
# shannon_genus
# 
# #El índice preterminado de la función de diversidad es el índice de Shannon, por lo que index = "shannon" se puede omitir
# shannon_genus<-diversity(abund_table, MARGIN=1)
# shannon_genus
# 
# #Dado que el cálculo del índice de Shannon es relativo a las filas. El argumento MARGIN=1 también se puede omitir
# shannon_genus<-diversity(abund_table)
# shannon_genus
# 
# #Ahora ilustraremos el cálculo por códigos del plan R utilizando la fórmula del índice de diversidad de Shannon-Wiener. Dado que en la fórmula anterior p_i es la proporción de individuos (o abundancia relativa) de la especie i en la comunidad, usamos la función decostand() en el paquete vegano para convertir los datos de conteo en cada muestra en proporciones.
# abund_table_total<-decostand(abund_table, MARGIN=1, method="total")
# 
# #Los aumentos de MARGIN = 1 significan "filas", MARGIN = 2 significan "columnas" de los datos del objeto tipo matriz (en este caso, abund_table); method = “total” significa dividir por el total del margen (MARGIN = 1 también es el valor predeterminado). Aplicando la función decostand se obtuvo la proporción de individuos en las muestras que pertenecen a cada gen. Luego podemos usar la fórmula del índice Shannon-Wiener para calcular el índice como se muestra a continuación.
# abund_table_p_lnp<-abund_table_total*log(abund_table_total)
# #sumamos los valores p_i*ln (p_i) y los multiplicamos por -1
# rowSums(abund_table_p_lnp, na.rm=TRUE)*-1
# 
# #En estos dos enfoques se obtienen los mismos resultados.
# 
# #----------------------------------------Índice de Simpson
# #Podemos usar la función diversity() en el paquete Vegan o la función R simple para calcular de índice de Simpson.
# 
# simp_genus<-diversity(abund_table, "simpson")
# simp_genus
# 
# #Los siguiente códigos de R utilizan la fórmula del índice de Simpson para calcularlo. Los datos de conteo deben convertirse en proporsiones antes de realizar el formulario.
# 
# #Usando decostand para convertir los datos en proporciones
# abund_table_total<-decostand(abund_table, MARGIN=1, method="total")
# #Elevamos al cuadrado
# abund_table_total_p2<-abund_table_total^2
# #Sumamos las filas
# 1-rowSums(abund_table_total_p2, na.rm = TRUE)
# 
# #El índice inverso de Simpson se puede calcular usando "invsimpson"
# 
# inv_simp<-diversity(abund_table, "invsimpson")
# inv_simp
# 
# #---------------------------------Índice de Equidad de Pielou
# 
# #El índice de Pielou se puede calcular utilizando las funciones specnumber() y diversity() basado en su fórmula que presentamos arriba
# S<-specnumber(abund_table)
# H<-diversity(abund_table, "shannon")
# J<-H/log(S)
# J
# 
# #--------------------------------Hacer un dataframe de índices de diversidad
# 
# #Podemos hacer un dataframe (marco de datos) para la cantidad de géneros usando los siguientes códigos R
# 
# #Hacer un marco de datos de número de géneros
# N<-specnumber(abund_table)
# df_N<-data.frame(sample=names(N), value=N, measure=rep("Number", length(N))) 
# 
# #El siguiente código se usa para crear un dataframe para el índice Chao1
# 
# #Hacer un marco de datos de número de géneros
# CH<-estimateR(abund_table)[2,]
# df_CH<-data.frame(sample=names(N), value=CH, measure=rep("Chao1", length(CH))) 
# 
# #El siguiente código se usa para crear un dataframe para el índice de Shannon
# 
# H<-diversity(abund_table, "shannon")
# df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon", length(H)))
# 
# #El siguiente código se usa para crear un dataframe para el índice de Simpson
# 
# df_simp<-data.frame(sample=names(simp_genus),value=simp_genus, measure=rep("Simpson",length(simp_genus)))
# 
# #El siguiente código se usa para crear un dataframe para el índice de Pielou
# 
# df_J<-data.frame(sample=names(J),value=J,measure=rep("Pielou", length(J)))
# 
# #Ahora podemos combinar todos los dataframe para usarlos después
# 
# df<-rbind(df_N, df_CH, df_H, df_simp, df_J)
# rownames(df)<-NULL
# df
# 
# 

