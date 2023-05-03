## ANALIZAMOS EL SEGUNDO CONJUNTO DE DATOS POR SEPARADO (3 CATEGORIAS POR TIPO DE CULTIVO)

library("phyloseq") 
library("stringi")
library("ggplot2")

# kraken-biom taxonomic_profiling/kraken_results/* --fmt json -o fresa_kraken_all.biom

## Hubo un problema con la creacion de BIOM,que no unia las nuevas muestras,lo cual se soluciono creando el BIOM de la siguiente manera

#~/GIT/Tesis_Maestria/Data/fresa_solena/Data_all/taxonomic_profiling/kraken_results
# $ conda activate metagenomics
# $ kraken-biom *txt --fmt json -o fresa_kraken.biom

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data2")
fresa_kraken <- import_biom("fresa_kraken.biom")
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
metadata_fresa <- read.csv2("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data2/metadata.csv",header =  TRUE, row.names = 1, sep = ",")
rownames(metadata_fresa) <- sample_names(fresa_kraken)
## unir los metadatos al objeto phyloseq
fresa_kraken@sam_data <- sample_data(metadata_fresa) 

## diversidad alfa
index = estimate_richness(fresa_kraken)
index
plot_richness(physeq = fresa_kraken, measures = c("Observed","Chao1","Shannon","simpson"))
plot_richness(physeq = fresa_kraken, measures = c("Observed","Chao1","Shannon","simpson"),x = "category", color = "category") 


## diversidad beta
percentages <- transform_sample_counts(fresa_kraken, function(x) x*100 / sum(x) )
head(percentages@otu_table@.Data)
meta_ord <- ordinate(physeq = percentages, method = "NMDS", distance = "bray") 
plot_ordination(physeq = percentages, ordination = meta_ord, color = "category") +
  geom_text(mapping = aes(label = colnames(fresa_kraken@otu_table@.Data)), size = 3, vjust = 1.5)

##### ERROR DATOS DE COVARIANZA

#####################################################################################################################

## StackBar

percentages_df <- psmelt(percentages)

# Ahora vamos a ordenar el data frame, para que nos quede en el orden que queremos graficar
percentages_df$Sample<-as.factor(percentages_df$Sample)
percentages_df$category<-as.factor(percentages_df$category)
# Ordenamos respecto a categoria 
percentages_df<-percentages_df[order(percentages_df$category,percentages_df$Sample),]

ggplot(data=percentages_df, aes_string(x='Sample', y='Abundance', fill='Phylum' ,color='category'))  +
  scale_colour_manual(values=c('cyan','pink','yellow')) +
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



