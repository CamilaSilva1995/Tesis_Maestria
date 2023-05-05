## SE ANALIZAN LA TOTALIDAD DELOS DATOS (85, 80 despues del filtro de calidad)

library("phyloseq") 
library("stringi")
library("ggplot2")

# kraken-biom taxonomic_profiling/kraken_results/* --fmt json -o fresa_kraken_all.biom

## Hubo un problema con la creacion de BIOM,que no unia las nuevas muestras,lo cual se soluciono creando el BIOM de la siguiente manera

#~/GIT/Tesis_Maestria/Data/fresa_solena/Data_all/taxonomic_profiling/kraken_results
# $ conda activate metagenomics
# $ kraken-biom *txt --fmt json -o fresa_kraken.biom

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data_all")
fresa_kraken <- import_biom("fresa_kraken_all.biom")
class(fresa_kraken)
## como tenemos diferencias de longitud en los nombres de las muestras, al cortar los nombres, quedan los datos nuevos con un punto al final, lo  cual para ciertos procedimientos a R no le gusta, por lo tanto a continuacion reemplazamos el punto por una X y asi no tener inconvenientes mas adelante 
## usando lalibreria stringi
sample_names(fresa_kraken)<-stri_replace_all_regex(sample_names(fresa_kraken),'\\.','A')
## renombrar las columnas del tax_table
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
## quitar los primeros caracteres de los nombres del tax_table
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data,4,100)
## recortar los nombres de las muestras
colnames(fresa_kraken@otu_table@.Data) <- substr(colnames(fresa_kraken@otu_table@.Data),1,6)
## cargar los metadatos
metadata_fresa <- read.csv2("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data_all/metadata.csv",header =  TRUE, row.names = 1, sep = ",")
rownames(metadata_fresa) <- sample_names(fresa_kraken)
## unir los metadatos al objeto phyloseq
fresa_kraken@sam_data <- sample_data(metadata_fresa) 

# filtro de calidad 
samples_to_remove <- c("MP2079","MP2080","MP2088","MP2109","MP2137")
samples<-!(sample_names(fresa_kraken) %in% samples_to_remove)
sample_names(fresa_kraken)
names(samples)<-as.character(sample_names(fresa_kraken))
samples
fresa_kraken_fil <- prune_samples(samples, fresa_kraken) #me esta quitando mas muestras de la que debe
nsamples(fresa_kraken) # 85
nsamples(fresa_kraken_fil) # 80
sample_names(fresa_kraken_fil) # 80#  

## diversidad alfa
index = estimate_richness(fresa_kraken_fil)
index
plot_richness(physeq = fresa_kraken_fil, measures = c("Observed","Chao1","Shannon","simpson"))
plot_richness(physeq = fresa_kraken_fil, measures = c("Observed","Chao1","Shannon","simpson"),x = "Category", color = "Category") 
ggsave("AllData_Alfa_diversidad.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

## diversidad beta
percentages_fil <- transform_sample_counts(fresa_kraken_fil, function(x) x*100 / sum(x) )
head(percentages_fil@otu_table@.Data)
meta_ord_fil <- ordinate(physeq = percentages_fil, method = "NMDS", distance = "bray") 
plot_ordination(physeq = percentages_fil, ordination = meta_ord_fil, color = "category") +
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)
ggsave("AllData_Beta_diversidad.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

#####################################################################################################################

## StackBar

percentages_fil <- transform_sample_counts(fresa_kraken_fil, function(x) x*100 / sum(x) )
percentages_df <- psmelt(percentages_fil)

# Ahora vamos a ordenar el data frame, para que nos quede en el orden que queremos graficar
percentages_df$Sample<-as.factor(percentages_df$Sample)
percentages_df$Category<-as.factor(percentages_df$Category)
# Ordenamos respecto a categoria 
percentages_df<-percentages_df[order(percentages_df$Category,percentages_df$Sample),]

##### Ahora vamos a crear una nueva variable, con los nombres ordenados de las muestras (New_sample)

crear_vector<-function(df,vec){
  t<-(as.integer(log10(length(vec))))+1 # numero de digitos del vector, se usa para calcular el numero de ceros 
  i<-1 # contador de muestra
  new_name<-c() # vector nombres nuevos
  new<-c() # 
  for (sample_i in unique(vec)){
    n <- count(df[vec==sample_i,]) #número de OTUs por muestra
    n <- as.integer(n)
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
  scale_colour_manual(values=c('white','black','cyan','pink','yellow')) +
  geom_bar(aes(), stat="identity", position="stack") +
  #scale_x_discrete(limits = rev(levels(percentages_df$Category))) +
  labs(title = "Abundance", x='Sample', y='Abundance', color = 'category') +
  theme(legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title=element_text(size=8, face = "bold"),
        legend.text=element_text(size=6),
        text = element_text(size=12),
        axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5))
ggsave("AllData_StackBar.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")




