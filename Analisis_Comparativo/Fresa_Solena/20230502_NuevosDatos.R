## ANALISIS PARA Data3... crop_vs_native_rawCounts_S
library("phyloseq") 
library("stringi")
library("ggplot2")

## Primero se debe crear el archivo BIOM
## kraken-biom taxonomic_profiling/kraken_results/* --fmt json -o crop_vs_native_kraken.biom

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data3")
fresa_kraken <- import_biom("crop_vs_native_kraken.biom")
class(fresa_kraken)
## como tenemos diferencias de longitud en los nombres de las muestras, al cortar los nombres, quedan los datos nuevos con un punto al final, lo  cual para ciertos procedimientos a R no le gusta, por lo tanto a continuacion reemplazamos el punto por una X y asi no tener inconvenientes mas adelante 
## usando la libreria stringi
sample_names(fresa_kraken)<-stri_replace_all_regex(sample_names(fresa_kraken),'\\.','A')
## renombrar las columnas del tax_table
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
## quitar los primeros caracteres de los nombres del tax_table
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data,4,100)

## cargar los metadatos
metadata_fresa <- read.csv2("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data3/metadata.csv",header =  TRUE, row.names = 1, sep = ",")
rownames(metadata_fresa) <- sample_names(fresa_kraken)

## unir los metadatos al objeto phyloseq
fresa_kraken@sam_data <- sample_data(metadata_fresa)

## Creamos una columna extra en sam_data por necesidad de funcionamiento de mas adelante
fresa_kraken@sam_data$Treatment<-c(rep('NA', nrow(metadata_fresa)))
colnames(fresa_kraken@sam_data)<-c('category','Treatment')
#View(fresa_kraken)

## Diversidad alfa
index = estimate_richness(fresa_kraken)
index
plot_richness(physeq = fresa_kraken, measures = c("Observed","Chao1","Shannon","simpson"))
plot_richness(physeq = fresa_kraken, measures = c("Observed","Chao1","Shannon","simpson"),x = "category", color = "category") 
ggsave("crop_vs_native_rawCounts_S_Alfa_Diversidad.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")


## Diversidad beta
percentages <- transform_sample_counts(fresa_kraken, function(x) x*100 / sum(x) )

meta_ord <- ordinate(physeq = percentages, method = "NMDS", distance = "bray") 
plot_ordination(physeq = percentages, ordination = meta_ord, color = "category") +
  geom_text(mapping = aes(label = colnames(fresa_kraken@otu_table@.Data)), size = 3, vjust = 1.5)
ggsave("crop_vs_native_rawCounts_S_Beta_Diversidad.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

## Barras de abundancia (StackBar)

## se convierte (percentages) a dataframe
percentages_df <- psmelt(percentages)
## Ahora vamos a ordenar el data frame, para que nos quede en el orden que queremos graficar
percentages_df$Sample<-as.factor(percentages_df$Sample)
percentages_df$category<-as.factor(percentages_df$category)
## Ordenamos respecto a la variable 'category' 
percentages_df<-percentages_df[order(percentages_df$category,percentages_df$Sample),]

## Funcion para crear vector de nueva variable con los nombres de las muestras ordenados
# crear_vector<-function(df,vec){
#   t<-(as.integer(log10(length(vec))))+1 # numero de digitos del vector, se usa para calcular el numero de ceros 
#   i<-1 # contador de muestra
#   new_name<-c() # vector nombres nuevos
#   new<-c() # 
#   for (sample_i in unique(vec)){
#     n <- count(df[vec==sample_i,]) #número de OTUs por muestra
#     n <- as.integer(n)
#     new <-c(rep(i, n)) # repetir i, n veces
#     m<-t-((as.integer(log10(i)))+1) # i es el numero de muestra, queremos "00i"
#     # acorde a la longitud de t es el numero de ceros para los nuevos nombres
#     # ejemplo t=3 queremos dos ceros 001
#     ceros<-paste((rep(0,m)),collapse = '')
#     new<-paste(paste('A',ceros,sep=""),new,sep = "") # nuevos nombres (por separado) "A0001"
#     # cat(i,as.integer(log10(i))+1,m,ceros,new,"\n")
#     new_name <- c(new_name,new) # vector acomula nuevos nombres 
#     i<-i+1
#   }
#   new_vec<-paste(new_name, vec, sep="_")
#   return(new_vec)
# }


#######################################################
df=percentages_df
vec=percentages_df$Sample
sample_i= '5QBDM2ESPPOOBXX0490'

  t<-(as.integer(log10(length(vec))))+1 # numero de digitos del vector, se usa para calcular el numero de ceros
  i<-1 # contador de muestra
  new_name<-c() # vector nombres nuevos
  new<-c() #
  #for (sample_i in unique(vec)){
    n <- count(df[vec==sample_i,]) #número de OTUs por muestra
    n <- as.integer(n)
    new <-c(rep(i, n)) # repetir i, n veces
    m<-t-((as.integer(log10(i)))+1) # i es el numero de muestra, queremos "00i"
    # acorde a la longitud de t es el numero de ceros para los nuevos nombres
    # ejemplo t=3 queremos dos ceros 001
    ceros<-paste((rep(0,m)),collapse = '')
    new<-paste(paste('A',ceros,sep=""),new,sep = "") # nuevos nombres (por separado) "A0001"
    # cat(i,as.integer(log10(i))+1,m,ceros,new,"\n")
    new_name <- c(new_name,new) # vector acomula nuevos nombres
    i<-i+1
  #}
#######################################################



## Ahora percentages_df$new_sample está ordenada alfabéticamente pero sigue el orden de las categorías
percentages_df$new_sample<-crear_vector(percentages_df,percentages_df$Sample)

## grafica barras de abundancia
ggplot(data=percentages_df, aes_string(x='Sample', y='Abundance', fill='Phylum' ,color='category'))  +
  scale_colour_manual(values=c('cyan','pink')) +
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

ggsave("crop_vs_native_rawCounts_S_StackBar.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")



