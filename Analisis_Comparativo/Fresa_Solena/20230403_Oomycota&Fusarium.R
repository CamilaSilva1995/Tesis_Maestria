## DIVERSIDADES ALFA Y BETA 
#Phylum - Oomycota
#Genus - Fusarium


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

## Aglomeramos Fusarium
glom_Fusarium <- subset_taxa(merge_Eukaryota, Genus == 'Fusarium')
glom_Fusarium_df <- psmelt(glom_Fusarium)

## se calcula diversidad alfa con el glom de fusarium 
plot_alpha_Fusarium <- plot_richness(physeq = glom_Fusarium, measures = c("Observed","Chao1","Shannon","simpson"),x = 'Treatment', color = 'Treatment') 
plot_alpha_Fusarium

## se calcula la beta diversidad para Fusarium
## sacamos los porcentajes
percentages_Fusarium <- transform_sample_counts(glom_Fusarium, function(x) x*100 / sum(x) )
percentages_Fusarium_df <- psmelt(percentages_Fusarium)

meta_ord_Fusarium <- ordinate(physeq = percentages_Fusarium, method = "NMDS", distance = 'bray') 
plot_beta_Fusarium <- plot_ordination(physeq = percentages_Fusarium, ordination = meta_ord_Fusarium, color = 'Treatment') +
  geom_text(mapping = aes(label = colnames(glom_Fusarium@otu_table@.Data)), size = 3, vjust = 1.5)
plot_beta_Fusarium

## graficamos barras de abundancia 
plot_barras_Fusarium <- ggplot(data = percentages_Fusarium_df, mapping = aes(y= Abundance, x = Sample, fill = Species , color = Treatment))+
  geom_bar(position = "stack", stat = "identity", size=1) +
  scale_color_manual(values=c('white','black')) +
  ylab(label = "Fusarium-Abundance") +
  theme(legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title=element_text(size=8, face = "bold"),
        legend.text=element_text(size=10),
        text = element_text(size=12),
        axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5))
plot_barras_Fusarium

## Aglomeramos Oomycota
glom_Oomycota <- subset_taxa(merge_Eukaryota, Phylum == 'Oomycota')
glom_Oomycota_df <- psmelt(glom_Oomycota)

## se calcula diversidad alfa con el glom de Oomycota 
plot_alpha_Oomycota <- plot_richness(physeq = glom_Oomycota, measures = c("Observed","Chao1","Shannon","simpson"),x = 'Treatment', color = 'Treatment') 
plot_alpha_Oomycota

## se calcula la beta diversidad para Oomycota
## sacamos los porcentajes
percentages_Oomycota <- transform_sample_counts(glom_Oomycota, function(x) x*100 / sum(x) )
percentages_Oomycota_df <- psmelt(percentages_Oomycota)

meta_ord_Oomycota <- ordinate(physeq = percentages_Oomycota, method = "NMDS", distance = 'bray') 
plot_beta_Oomycota <- plot_ordination(physeq = percentages_Oomycota, ordination = meta_ord_Oomycota, color = 'Treatment') +
  geom_text(mapping = aes(label = colnames(glom_Oomycota@otu_table@.Data)), size = 3, vjust = 1.5)
plot_beta_Oomycota

## graficamos barras de abundancia 
plot_barras_Oomycota <- ggplot(data = percentages_Oomycota_df, mapping = aes(y= Abundance, x = Sample, fill = Species , color = Treatment))+
  geom_bar(position = "stack", stat = "identity", size=1) +
  scale_color_manual(values=c('white','black')) +
  ylab(label = "Oomycota-Abundance") +
  theme(legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title=element_text(size=8, face = "bold"),
        legend.text=element_text(size=10),
        text = element_text(size=12),
        axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5))
plot_barras_Oomycota




