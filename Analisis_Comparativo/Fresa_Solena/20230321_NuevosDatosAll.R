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


## diversidad beta
percentages_fil <- transform_sample_counts(fresa_kraken_fil, function(x) x*100 / sum(x) )
head(percentages_fil@otu_table@.Data)
meta_ord_fil <- ordinate(physeq = percentages_fil, method = "NMDS", distance = "bray") 
plot_ordination(physeq = percentages_fil, ordination = meta_ord_fil, color = "Treatment") +
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)

plot_ordination(physeq = percentages_fil, ordination = meta_ord_fil, color = "Category") +
  geom_text(mapping = aes(label = colnames(fresa_kraken_fil@otu_table@.Data)), size = 3, vjust = 1.5)

#####################################################################################################################

