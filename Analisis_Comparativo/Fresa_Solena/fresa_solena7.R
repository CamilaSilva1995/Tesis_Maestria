
library("phyloseq")
library("ggplot2")
library("vegan")
library("RColorBrewer")
library("stringi")

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena")
outpath = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"

### Cargado de datos originales 
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


# merge_Eukaryota<-subset_taxa(fresa_kraken_fil,Genus=="Fusarium")
# subconjunto cortado a nivel de genero, objeto phyloseq
Genus_glom <- tax_glom(fresa_kraken_fil, taxrank = 'Genus')
# se pasa a dataframe
Genus_glom_df = Genus_glom@otu_table
# se guarda como csv
write.csv(Genus_glom_df, "/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/fresa_genus.csv")
#write_tsv(Genus_glom_df, "/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/fresa_genus.csv")
write.table(Genus_glom_df, file="/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/fresa_genus.tsv", quote=FALSE, sep='\t')

# 















fresa_genus <- import_biom("fresa_genus.biom")



