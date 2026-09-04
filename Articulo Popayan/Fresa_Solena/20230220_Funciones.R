## CREAR UNA FUNCION QUE GRAFIQUE TODAS LAS OPCIONES

library("phyloseq")
library("ggplot2")
library("vegan")
library("RColorBrewer")

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
percentages_df <- psmelt(percentages_fil)

## Creamos dos funciones con las cuales automatizamos la creacion de los subconjuntos por nivel taxonomico y la creacion de las graficas de beta diversidad y abundacias de barras por cada subconjunto 

# input(objeto phyloseq, nivel taxonomico, 'Treatment')
## objeto phyloseq = fresa_kraken_fil
## crear subset 
## sacar los porcentajes 
## plot de abundancias con barras
## plot de beta diversidad
# output(plot1,plot2)

## objeto total
## objeto de bacterias
## objeto de hongos
## Con glom aglomeras por el rango que quieres, y ahora si haces los porcentajes y las graficas

## Subconjunto de "Eukaryota"
merge_Eukaryota<-subset_taxa(fresa_kraken_fil,Kingdom=="Eukaryota")
## Subconjunto de "Bacteria"
merge_Bacteria<-subset_taxa(fresa_kraken_fil,Kingdom=="Bacteria")

## inputs: phyloseq_objet ,character:taxrank,character: Atributo 
glomToGraph<-function(phy,tax,attribute,abundance_percentage){
  ## creamos el subconjunto dependiendo del linaje taxonomico deseado
  glom <- tax_glom(phy, taxrank = tax)
  ## sacamos los porcentajes
  percentages <- transform_sample_counts(glom, function(x) x*100 / sum(x) )
  percentages_df <- psmelt(percentages)
  ## Grafica de diversidad en barras
  plot_barras <- ggplot(data=percentages_df, aes_string(x='Sample', y='Abundance', fill=tax ,color = attribute))+
    geom_bar(aes(), stat="identity", position="stack") +
    theme(text = element_text(size = 5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ## Ahora queremos un glom por porcentajes
  percentages_df$tax<-percentages_df[,ncol(percentages_df)]
  percentages_df$tax[percentages_df$Abundance < abundance_percentage] <- "abundance_percentage" 
  percentages_df$tax <- as.factor(percentages_df$tax)
  plot_percentage <- ggplot(data=percentages_df, aes_string(x='Sample', y='Abundance', fill='tax' ,color = attribute))+
    geom_bar(aes(), stat="identity", position="stack") +
    theme(text = element_text(size = 5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  return(list(plot_barras,plot_percentage))
}

## Esta funcion recibe un objeto phyloseq, un rango taxonomico, una variable de atributo y una distancia de beta diversidad, 
Beta_diversity <- function(phy,tax,attribute,distance){
  ## creamos el subconjunto dependiendo del linaje taxonomico deseado
  merge<-subset_taxa(phy, taxrank = tax)
  ## sacamos los porcentajes
  percentages <- transform_sample_counts(merge, function(x) x*100 / sum(x) )
  percentages_df <- psmelt(percentages)
  # head(percentages_df)
  # return(head(percentages_df))
  ## Beta diversidad 
  meta_ord <- ordinate(physeq = percentages, method = "NMDS", distance = distance) 
  plot_beta <- plot_ordination(physeq = percentages, ordination = meta_ord, color = attribute) +
    geom_text(mapping = aes(label = colnames(phy@otu_table@.Data)), size = 3, vjust = 1.5)
  return(plot_beta)
}

Alpha_diversity <- function(phy,tax,attribute,distance){
  ## creamos el subconjunto dependiendo del linaje taxonomico deseado
  merge<-subset_taxa(phy, taxrank = tax)
  ## sacamos los porcentajes
  percentages <- transform_sample_counts(merge, function(x) x*100 / sum(x) )
  percentages_df <- psmelt(percentages)
  # head(percentages_df)
  # return(head(percentages_df))
  ## Beta diversidad 
  meta_ord <- ordinate(physeq = percentages, method = "NMDS", distance = distance) 
  plot_alpha <- plot_ordination(physeq = percentages, ordination = meta_ord, color = attribute) +
    geom_text(mapping = aes(label = colnames(phy@otu_table@.Data)), size = 3, vjust = 1.5)
  return(plot_alpha)
}


Beta_Phylum <- Beta_diversity(fresa_kraken_fil,'Phylum','Treatment',"bray")
ggsave("Beta_Phylum.png", plot = Beta_Phylum, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
Barras_Phylum <- glomToGraph(merge_Eukaryota,'Phylum' , 'Treatment', 10.0)
Barras_Phylum[1] # normal
ggsave("Barras_Phylum_Eukaryota.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
Barras_Phylum[2] # con porcentaje de abundancia
ggsave("Barras_Phylum_Eukaryota_percentageAbundance.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
Barras_Phylum <- glomToGraph(merge_Bacteria,'Phylum' , 'Treatment', 10.0)
Barras_Phylum[1] # normal
ggsave("Barras_Phylum_Bacteria.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
Barras_Phylum[2] # con porcentaje de abundancia
ggsave("Barras_Phylum_Bacteria_percentageAbundance.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Beta_Family <- Beta_diversity(fresa_kraken_fil,'Family','Treatment',"bray")
ggsave("Beta_Family.png", plot = Beta_Family, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
Barras_Family <- glomToGraph(merge_Eukaryota,'Family' , 'Treatment')
ggsave("Barras_Family_Eukaryota.png", plot = Barras_Family, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
Barras_Family <- glomToGraph(merge_Bacteria,'Family' , 'Treatment')
ggsave("Barras_Family_Bacteria.png", plot = Barras_Family, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Beta_Genus <- Beta_diversity(fresa_kraken_fil,'Genus','Treatment',"bray")
ggsave("Beta_Genus.png", plot = Beta_Genus, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
Barras_Genus <- glomToGraph(merge_Eukaryota,'Genus' , 'Treatment')
ggsave("Barras_Genus_Eukaryota.png", plot = Barras_Genus, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 20, height = 10, dpi = 300, units = "cm")
Barras_Genus <- glomToGraph(merge_Bacteria,'Genus' , 'Treatment')
ggsave("Barras_Genus_Bacteria.png", plot = Barras_Genus, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 10, height = 5, dpi = 300, units = "cm")

Beta_Species <- Beta_diversity(fresa_kraken_fil,'Species','Treatment',"bray")
ggsave("Beta_Species.png", plot = Beta_Species, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
Barras_Species <- glomToGraph(merge_Eukaryota,'Species' , 'Treatment')
ggsave("Barras_Species_Eukaryota.png", plot = Barras_Species, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 10, height = 5, dpi = 300, units = "cm")
Barras_Species <- glomToGraph(merge_Bacteria,'Species' , 'Treatment')
ggsave("Barras_Species_Bacteria.png", plot = Barras_Species, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 10, height = 5, dpi = 300, units = "cm")








