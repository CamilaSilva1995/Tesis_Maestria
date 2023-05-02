## QUEREMOS GRAFICAR LOS 80 DATOS ORDENADOS POR GRUPOS

library("phyloseq")
library("ggplot2")
library("vegan")
library("RColorBrewer")
library("stringi")
library("forcats")
library("dplyr")

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data_all")
fresa_kraken <- import_biom("fresa_kraken_all.biom")
class(fresa_kraken)
## como tenemos diferencias de longitud en los nombres de las muestras, al cortar los nombres, quedan los datos nuevos con un punto al final, lo  cual para ciertos procedimientos a R no le gusta, por lo tanto a continuacion reemplazamos el punto por una X y asi no tener inconvenientes mas adelante
## usando la libreria stringi
sample_names(fresa_kraken)<-stri_replace_all_regex(sample_names(fresa_kraken),'\\.kraken2.report','')
## renombrar las columnas del tax_table
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
## quitar los primeros caracteres de los nombres del tax_table
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data,4,100)
## cargar los metadatos
metadata_fresa <- read.csv2("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data_all/metadata.csv",header =  TRUE, row.names = 1, sep = ",")

#metadata_fresa
## unir los metadatos al objeto phyloseq
fresa_kraken@sam_data <- sample_data(metadata_fresa)

# filtro de calidad
samples_to_remove <- c("MP2079","MP2080","MP2088","MP2109","MP2137")
samples<-!(sample_names(fresa_kraken) %in% samples_to_remove)
names(samples)<-as.character(sample_names(fresa_kraken))
fresa_kraken_fil <- prune_samples(samples, fresa_kraken)
nsamples(fresa_kraken) # 85
nsamples(fresa_kraken_fil) # 80

percentages_fil <- transform_sample_counts(fresa_kraken_fil, function(x) x*100 / sum(x) )
percentages_df <- psmelt(percentages_fil)

percentages_df$Sample<-as.factor(percentages_df$Sample)
percentages_df$Category<-as.factor(percentages_df$Category)
percentages_df<-percentages_df[order(percentages_df$Category,percentages_df$Sample),]

# percentages_df<-percentages_df[order(percentages_df$Category,percentages_df$Sample),]
# #x$name <- factor(x$name, levels = x$name[order(x$val)])
# #percentages_df$Sample<-factor(percentages_df$Sample, levels = percentages_df$Sample[order(percentages_df$Category)])
#
#
# #percentages_df<-percentages_df %>%    # First sort by val. This sort the dataframe but NOT the factor levels
# #  mutate(Sample=factor(Sample, levels=Sample)) # This trick update the factor levels
# #head(percentages_df)
#
# #x$name <- factor(x$name, levels = x$name[order(x$val)])
# #percentages_df$Sample <- factor(percentages_df$Sample, levels =percentages_df$Sample[order(percentages_df$Sample)])
#
# percentages_df %>% mutate(Category = factor(Category, levels = c("cropHealthy","cropWilted","nearCrop", "forestDegraded","uncultivatedLand")))



crear_vector<-function(df,vec){
  t<-(as.integer(log10(length(vec))))+1
  i<-1
  new_name<-c()
  new<-c()
  for (sample_i in unique(vec)){
    n <- count(df[vec==sample_i,]) #número de OTUs por muestra
    n <- as.integer(n)
    return(n) # YA NO FUNCIONA EL CONTEO DE OTUS
    new <-c(rep(i, n))
    m<-t-((as.integer(log10(i)))+1)
    ceros<-paste((rep(0,m)),collapse = '')
    new<-paste(paste('A',ceros,sep=""),new,sep = "")
    # cat(i,as.integer(log10(i))+1,m,ceros,new,"\n")
    new_name <- c(new_name,new)
    i<-i+1
  }
  new_vec<-paste(new_name, vec, sep="_")
  return(new_vec)
}


percentages_df$new_sample<-crear_vector(percentages_df,percentages_df$Sample)

ggplot(data=percentages_df, aes_string(x='Sample', y='Abundance', fill='Phylum' ,color='Category'))  +
  scale_colour_manual(values=c('white','black','cyan','pink','yellow')) +
  geom_bar(aes(), stat="identity", position="stack") +
  #scale_x_discrete(limits = rev(levels(percentages_df$Category))) +
  labs(title = "Abundance", x='Sample', y='Abundance', color = 'Category') +
  theme(legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title=element_text(size=8, face = "bold"),
        legend.text=element_text(size=6),
        text = element_text(size=12),
        axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5))



#queremos crear una funcion que haga un nuevo vector de nombres con orden alfabetico
# new_name<-c()
# new<-c()
# for (i in 1:length(unique(little$Sample))){
#   n <- count(little[little$Sample==little$Sample[i],])
#   n <- as.integer(n)
#   for(j in 1:n){
#     new<-c(rep(i, n))
#     new<-paste('A0',new,sep="")
#     # print(n)
#     # print(new)
#   }
#   new_name <- c(new_name,new)
#   # print(new_name)
# }
# little$new_sample<-paste(new_name, little$Sample, sep="_")
# little

# 
###EJEMPLO PEQUEÑO
# 
percentages_df2<-read.csv2("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/percentages_df.csv",header =  TRUE, sep = ",")
head(percentages_df2)
#La función rbind en R, abreviatura de row-bind , se puede usar para combinar vectores, matrices y marcos de datos por filas.
little<-rbind(percentages_df2[1:7,],percentages_df2[277:280,],percentages_df2[237:244,],percentages_df2[213:216,])
unique(little$Sample)
little

#MP2047 debe ser blanco
#MP2048 debe ser blanco
#P3029 debe ser cyan
#P3358 debe ser cyan poco verde (7.13)
#P3066 debe ser negro de verde (9.2)  
#P3399 debe ser negro                 
# Salen ordenados de acuerdo a sample alfabeticamente
ggplot(data=little, aes_string(x=('Sample'), y='Abundance', fill='Phylum', colour='Category' ))+
  scale_colour_manual(values=c('white','black','cyan','pink','yellow')) +
  geom_bar(aes(), stat="identity", position="stack")

little<-little[order(little$Category,little$Sample),]



## QUEREMOS GRAFICAR LOS 80 DATOS ORDENADOS POR GRUPOS

library("phyloseq")
library("ggplot2")
library("vegan")
library("RColorBrewer")
library("stringi")
library("forcats")
library("dplyr")

setwd("/home/nelly/Tesis_Maestria/Data/fresa_solena/Data_all")
fresa_kraken <- import_biom("fresa_kraken_all.biom")
class(fresa_kraken)
## como tenemos diferencias de longitud en los nombres de las muestras, al cortar los nombres, quedan los datos nuevos con un punto al final, lo  cual para ciertos procedimientos a R no le gusta, por lo tanto a continuacion reemplazamos el punto por una X y asi no tener inconvenientes mas adelante
## usando la libreria stringi
sample_names(fresa_kraken)<-stri_replace_all_regex(sample_names(fresa_kraken),'\\.kraken2.report','')
## renombrar las columnas del tax_table
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
## quitar los primeros caracteres de los nombres del tax_table
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data,4,100)
## cargar los metadatos
metadata_fresa <- read.csv2("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data_all/metadata.csv",header =  TRUE, row.names = 1, sep = ",")

#metadata_fresa
## unir los metadatos al objeto phyloseq
fresa_kraken@sam_data <- sample_data(metadata_fresa)

# filtro de calidad
samples_to_remove <- c("MP2079","MP2080","MP2088","MP2109","MP2137")
samples<-!(sample_names(fresa_kraken) %in% samples_to_remove)
names(samples)<-as.character(sample_names(fresa_kraken))
fresa_kraken_fil <- prune_samples(samples, fresa_kraken)
nsamples(fresa_kraken) # 85
nsamples(fresa_kraken_fil) # 80

percentages_fil <- transform_sample_counts(fresa_kraken_fil, function(x) x*100 / sum(x) )
percentages_df <- psmelt(percentages_fil)

# Ahora vamos a ordenar el data frame, para que nos quede en el orden que queremos graficar
percentages_df$Sample<-as.factor(percentages_df$Sample)
percentages_df$Category<-as.factor(percentages_df$Category)
# Ordenamos respecto a categoria 
percentages_df<-percentages_df[order(percentages_df$Category,percentages_df$Sample),]

##### Ahora vamos a crear una nueva variable, con los nombres ordenados de las muestras (New_sample)
## Ejemplo este es un data frame input
#       OTU1 MD2056 
#       OTU2 MD2056     
#       OTU3 MD2056 
#       OTU1 MD2056 
#       OTU2 MD2066  
#       OTU3 MD2066
#       OTU1 MD2058 
#       OTU2 MD2058

# QUEREMOS este vector reordenado de salida
#       OTU1 A001_MD2056 
#       OTU2 A001_MD2056     
#       OTU3 A001_MD2056 
#       OTU1 A003_MD2056 
#       OTU2 A003_MD2066  
#       OTU3 A003_MD2066
#       OTU1 A002_MD2058 
#       OTU2 A002_MD2058


crear_vector<-function(df,vec){
  t<-(as.integer(log10(length(vec))))+1 # numero de digitos del vector, se usa para calcular el numero de ceros 
  i<-1 # contador de muestra
  new_name<-c() # vector nombres nuevos
  new<-c() # 
  for (sample_i in unique(vec)){
# sample_i<-"P3367"
   n <- count(df[vec==sample_i,]) #número de OTUs por muestra
    n <- as.integer(n)
    print(n)
    new <-c(rep(i, n)) # repetir i, n veces
    m<-t-((as.integer(log10(i)))+1) # i es el numero de muestra, queremos "00i"
                                    #acorde a la longitud de t es el numero de ceros para los nuevos nombres
                                    # ejemplo t=3 queremos dos ceros 001
    ceros<-paste((rep(0,m)),collapse = '')
    new<-paste(paste('A',ceros,sep=""),new,sep = "") # nuevos nombres (por separado) "A0001"
    # cat(i,as.integer(log10(i))+1,m,ceros,new,"\n")
    new_name <- c(new_name,new) # vector acomula nuevos nombres 
    i<-i+1
   }
  new_vec<-paste(new_name, vec, sep="_")
  return(new_vec)
}

### Ahora percentages_df$new_sample está ordenada alfabéticamente pero sigue el orden de las categorías
percentages_df$new_sample<-crear_vector(percentages_df,percentages_df$Sample)

ggplot(data=percentages_df, aes_string(x='new_sample', y='Abundance', fill='Phylum' ,color='Category'))  +
 # scale_colour_manual(values=c('white','black','cyan','pink','yellow')) +
  geom_bar(aes(), stat="identity", position="stack") +
  #scale_x_discrete(limits = rev(levels(percentages_df$Category))) +
  labs(title = "Abundance", x='Sample', y='Abundance', color = 'Category') +
  theme(legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title=element_text(size=8, face = "bold"),
        legend.text=element_text(size=6),
        text = element_text(size=12),
        axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5))



###EJEMPLO PEQUEÑO
# 
percentages_df2<-read.csv2("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/percentages_df.csv",header =  TRUE, sep = ",")
head(percentages_df2)
#La función rbind en R, abreviatura de row-bind , se puede usar para combinar vectores, matrices y marcos de datos por filas.
little<-rbind(percentages_df2[1:7,],percentages_df2[277:280,],percentages_df2[237:244,],percentages_df2[213:216,])
unique(little$Sample)
little

#MP2047 debe ser blanco
#MP2048 debe ser blanco
#P3029 debe ser cyan
#P3358 debe ser cyan poco verde (7.13)
#P3066 debe ser negro de verde (9.2)  
#P3399 debe ser negro                 
# Salen ordenados de acuerdo a sample alfabeticamente
ggplot(data=little, aes_string(x=('Sample'), y='Abundance', fill='Phylum', colour='Category' ))+
  scale_colour_manual(values=c('white','black','cyan','pink','yellow')) +
  geom_bar(aes(), stat="identity", position="stack")

little<-little[order(little$Category,little$Sample),]









