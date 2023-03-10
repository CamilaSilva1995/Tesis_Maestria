## DIVERSIDADES ALFA Y BETA FUSARIUM

library("phyloseq")
library("ggplot2")
library("vegan")
#library("BiodiversityR") 
library("RColorBrewer")
library("stringi")

#https://ciespinosa.github.io/AlphaDiversidad/rarefaccion.html

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
percentages_fil <- transform_sample_counts(fresa_kraken_fil, function(x) x*100 / sum(x) )
percentages_df <- psmelt(percentages_fil)

## Subconjunto de "Eukaryota"
merge_Eukaryota<-subset_taxa(fresa_kraken_fil,Kingdom=="Eukaryota")
## Subconjunto de "Bacteria"
merge_Bacteria<-subset_taxa(fresa_kraken_fil,Kingdom=="Bacteria")

## Aglomeramos Fusarium
glom <- tax_glom(merge_Eukaryota, taxrank = 'Genus')
glom_Fusarium <- subset_taxa(glom, Genus == 'Fusarium')

## se calcula diversidad alfa con el glom de fusarium 
plot_alpha_Fusarium <- plot_richness(physeq = glom_Fusarium, measures = c("Observed","Chao1","Shannon","simpson"),x = 'Treatment', color = 'Treatment') 
plot_alpha_Fusarium

## se calcula la beta diversidad para Fusarium
## sacamos los porcentajes
percentages_Fusarium <- transform_sample_counts(glom_Fusarium, function(x) x*100 / sum(x) )
percentages_Fusarium_df <- psmelt(percentages_Fusarium)
meta_ord_Fusarium <- ordinate(physeq = percentages_Fusarium, method = "NMDS", distance = 'bray') 
plot_bet_Fusariuma <- plot_ordination(physeq = percentages_Fusarium, ordination = meta_ord_Fusarium, color = 'Treatment') +
  geom_text(mapping = aes(label = colnames(phy@otu_table@.Data)), size = 3, vjust = 1.5)
plot_beta_Fusarium




## Aglomeramos Actinobacteria
#Actinobacteria
glom <- tax_glom(merge_Bacteria, taxrank = 'Phylum')
glom_Actinobacteria <- subset_taxa(glom, Phylum == 'Actinobacteria')

## se calcula diversidad alfa con el glom de Actinobacteria 
plot_alpha_Actinobacteria <- plot_richness(physeq = glom_Actinobacteria, measures = c("Observed","Chao1","Shannon","simpson"),x = 'Treatment', color = 'Treatment') 
plot_alpha_Actinobacteria

## se calcula la beta diversidad para Actinobacteria
## sacamos los porcentajes
percentages_Actinobacteria <- transform_sample_counts(glom_Actinobacteria, function(x) x*100 / sum(x) )
percentages_Actinobacteria_df <- psmelt(percentages_Actinobacteria)
meta_ord_Actinobacteria <- ordinate(physeq = percentages_Actinobacteria, method = "NMDS", distance = 'bray') 
plot_bet_Actinobacteriaa <- plot_ordination(physeq = percentages_Actinobacteria, ordination = meta_ord_Actinobacteria, color = 'Treatment') +
  geom_text(mapping = aes(label = colnames(phy@otu_table@.Data)), size = 3, vjust = 1.5)
plot_beta_Actinobacteria


