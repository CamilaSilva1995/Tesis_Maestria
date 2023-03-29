library("phyloseq")
library("ggplot2")
#library("vegan")
library("RColorBrewer")
library("stringi")
library(forcats)

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data_all")
fresa_kraken <- import_biom("fresa_kraken_all.biom")
class(fresa_kraken)
## como tenemos diferencias de longitud en los nombres de las muestras, al cortar los nombres, quedan los datos nuevos con un punto al final, lo  cual para ciertos procedimientos a R no le gusta, por lo tanto a continuacion reemplazamos el punto por una X y asi no tener inconvenientes mas adelante 
## usando lalibreria stringi
sample_names(fresa_kraken)<-stri_replace_all_regex(sample_names(fresa_kraken),'\\.kraken2.report','')
## renombrar las columnas del tax_table
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
## quitar los primeros caracteres de los nombres del tax_table
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data,4,100)
## recortar los nombres de las muestras
#colnames(fresa_kraken@otu_table@.Data) <- substr(colnames(fresa_kraken@otu_table@.Data),1,6)
## cargar los metadatos
metadata_fresa <- read.csv2("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data_all/metadata.csv",header =  TRUE, row.names = 1, sep = ",")
#rownames(metadata_fresa) <- sample_names(fresa_kraken)
#rownames(metadata_fresa) <- metadata_fresa$

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

merge_Eukaryota<-subset_taxa(fresa_kraken_fil,Kingdom=="Eukaryota")

glomToGraph<-function(phy,tax){
  ## creamos el subconjunto dependiendo del linaje taxonomico deseado
  glom <- tax_glom(phy, taxrank = tax)
  ## sacamos los porcentajes
  percentages <- transform_sample_counts(glom, function(x) x*100 / sum(x) )
  percentages_df <- psmelt(percentages)
  #return(nsamples(percentages_df))
  return(list(glom,percentages,percentages_df))
}

datos<-glomToGraph(merge_Eukaryota,'Phylum')
percentages_df<-datos[[3]]
length(percentages_df$Sample)


#percentages_df$Sample<-as.factor(percentages_df$Sample)
#percentages_df$Category<-as.factor(percentages_df$Category)

percentages_df<-percentages_df[order(percentages_df$Category,percentages_df$Sample),]
#x$name <- factor(x$name, levels = x$name[order(x$val)])
#percentages_df$Sample<-factor(percentages_df$Sample, levels = percentages_df$Sample[order(percentages_df$Category)])


#percentages_df<-percentages_df %>%    # First sort by val. This sort the dataframe but NOT the factor levels
#  mutate(Sample=factor(Sample, levels=Sample)) # This trick update the factor levels
#head(percentages_df)

#x$name <- factor(x$name, levels = x$name[order(x$val)])
#percentages_df$Sample <- factor(percentages_df$Sample, levels =percentages_df$Sample[order(percentages_df$Sample)])


ggplot(data=percentages_df, aes_string(x='Sample', y='Abundance', fill='Phylum' ,color='Category'))  +
  scale_colour_manual(values=c('white','black','cyan','pink','yellow')) +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_x_discrete(limits = rev(levels(percentages_df$Category))) +
  labs(title = "Abundance", x='Sample', y='Abundance', color = 'Category') +
  theme(legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title=element_text(size=8, face = "bold"),
        legend.text=element_text(size=6),
        text = element_text(size=12),
        axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5))






















