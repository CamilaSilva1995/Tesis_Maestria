## CREAR UNA FUNCION QUE GRAFIQUE TODAS LAS OPCIONES

library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("stringi")

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
outpath = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"

### Cargado de datos originales 
fresa_kraken <- import_biom("fresa_kraken.biom")
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data,4,100)
colnames(fresa_kraken@otu_table@.Data) <- substr(colnames(fresa_kraken@otu_table@.Data),1,6)
metadata_fresa <- read.csv2("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1/metadata.csv",header =  FALSE, row.names = 1, sep = ",")
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


##########Funciones 
##Crea los subconjuntos de datos 
## input
##       phy phyliseq total
##       tax rango al que queremos recortar
## output

glomToGraph<-function(phy,tax){
  ## creamos el subconjunto dependiendo del linaje taxonomico deseado
  glom <- tax_glom(phy, taxrank = tax)
  ## sacamos los porcentajes
  percentages <- transform_sample_counts(glom, function(x) x*100 / sum(x) )
  percentages_df <- psmelt(percentages)
  return(list(glom,percentages,percentages_df))
}
#-------------------------------------------------------------------------------
## Graficar abundancias stackbar

# input entra el percentages_df 
Abundance_barras <- function(phy,tax,attribute,abundance_percentage){
  ##llamar funcion de datos 
  Data <- glomToGraph(phy,tax)
  glom <- Data[[1]] #phyloseq
  percentages <- Data[[2]] #phyloseq
  percentages_df <- Data[[3]] # dataframe
  ## Graficamos para cada subconjunto las barras de abundancia
  plot_barras <- ggplot(data=percentages_df, aes_string(x='Sample', y='Abundance', fill=tax ,color=attribute)) + 
    scale_colour_manual(values=c('white','black')) +
    geom_bar(aes(), stat="identity", position="stack") +
    labs(title = "Abundance", x='Sample', y='Abundance', color = tax) +
    theme(legend.key.size = unit(0.2, "cm"),
          legend.key.width = unit(0.25,"cm"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title=element_text(size=8, face = "bold"),
          legend.text=element_text(size=6),
          text = element_text(size=12),
          axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5))
  percentages_df$tax<-percentages_df[,ncol(percentages_df)]
  percentages_df$tax[percentages_df$Abundance < abundance_percentage] <- "Others"
  percentages_df$tax <- as.factor(percentages_df$tax)
  plot_percentages <- ggplot(data=percentages_df, aes_string(x='Sample', y='Abundance', fill='tax' ,color=attribute))+
    scale_colour_manual(values=c('white','black')) +
    geom_bar(aes(), stat="identity", position="stack") +
    labs(title = "Abundance", x='Sample', y='Abundance', color = tax) +
    theme(legend.key.size = unit(0.3, "cm"),
          legend.key.width = unit(0.5,"cm"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title=element_text(size=10, face = "bold"),
          legend.text=element_text(size=8),
          text = element_text(size=12),
          axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5))
  return(list(plot_barras,plot_percentages))
}

#-----------------------------------------
#Graficar Beta Diversity
Beta_diversity <- function(phy,tax,attribute,distance){
  Data <- glomToGraph(phy,tax)
  glom <- Data[[1]]
  #CREAR UN GLOM AL 10%
  percentages <- Data[[2]]
  percentages_df <- Data[[3]]
  ## Beta diversidad 
  meta_ord <- ordinate(physeq = percentages, method = "NMDS", distance = distance) 
  plot_beta <- plot_ordination(physeq = percentages, ordination = meta_ord, color = attribute) +
    geom_text(mapping = aes(label = colnames(phy@otu_table@.Data)), size = 3, vjust = 1.5)
  return(plot_beta)
}

#Grafica Alfa Diversidad
Alpha_diversity <- function(phy,tax,attribute){
  ## llamamos la funcion que crea los dataset
  Data <- glomToGraph(phy,tax)
  glom <- Data[[1]]
  
  percentages <- Data[[2]]
  percentages_df <- Data[[3]]
  ## Alfa diversidad
  plot_alpha <- plot_richness(physeq = glom, measures = c("Observed","Chao1","Shannon","simpson"),x = attribute, color = attribute) 
  return(plot_alpha)
}


################################################################################
#### Guarda Graficas
################################################################################

############Eukarya#############################################################

#-----------Eukarya by Phylum 
Barras_Eukarya_Phylum <- Abundance_barras(merge_Eukaryota,'Phylum' , 'Treatment', 10.0)
Barras_Eukarya_Phylum[[1]] # normal
Barras_Eukarya_Phylum[[2]] # 10%
ggsave("Barras_Eukarya_Phylum.svg", plot = Barras_Eukarya_Phylum[[1]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
ggsave("Barras_Eukarya_Phylum_10.svg", plot = Barras_Eukarya_Phylum[[2]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Beta_Eukarya_Phylum<-Beta_diversity(merge_Eukaryota , 'Phylum' , 'Treatment', 'bray')
ggsave("Beta_Eukarya_Phylum.svg", plot = Beta_Eukarya_Phylum, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Alpha_Eukarya_Phylum<-Alpha_diversity(merge_Eukaryota , 'Phylum' , 'Treatment')
ggsave("Alpha_Eukarya_Phylum.svg", plot = Alpha_Eukarya_Phylum, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

#-----------Eukarya by Family 
Barras_Eukarya_Family <- Abundance_barras(merge_Eukaryota,'Family' , 'Treatment', 10.0)
Barras_Eukarya_Family[[1]] # normal
Barras_Eukarya_Family[[2]] # 10%
ggsave("Barras_Eukarya_Family.svg", plot = Barras_Eukarya_Family[[1]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
ggsave("Barras_Eukarya_Family_10.svg", plot = Barras_Eukarya_Family[[2]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Beta_Eukarya_Family<-Beta_diversity(merge_Eukaryota , 'Family' , 'Treatment', 'bray')
ggsave("Beta_Eukarya_Family.svg", plot = Beta_Eukarya_Family, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Alpha_Eukarya_Family<-Alpha_diversity(merge_Eukaryota , 'Family' , 'Treatment')
ggsave("Alpha_Eukarya_Family.svg", plot = Alpha_Eukarya_Family, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

#-----------Eukarya by Genus
Barras_Eukarya_Genus <- Abundance_barras(merge_Eukaryota,'Genus' , 'Treatment', 10.0)
Barras_Eukarya_Genus[[1]] # normal
Barras_Eukarya_Genus[[2]] # 10%
ggsave("Barras_Eukarya_Genus.svg", plot = Barras_Eukarya_Genus[[1]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
ggsave("Barras_Eukarya_Genus_10.svg", plot = Barras_Eukarya_Genus[[2]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Beta_Eukarya_Genus<-Beta_diversity(merge_Eukaryota , 'Genus' , 'Treatment', 'bray')
ggsave("Beta_Eukarya_Genus.svg", plot = Beta_Eukarya_Genus, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Alpha_Eukarya_Genus<-Alpha_diversity(merge_Eukaryota , 'Genus' , 'Treatment')
ggsave("Alpha_Eukarya_Genus.svg", plot = Alpha_Eukarya_Genus, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

#-----------Eukarya by Species 

Barras_Eukarya_Species <- Abundance_barras(merge_Eukaryota,'Species' , 'Treatment', 10.0)
Barras_Eukarya_Species[[1]] 
Barras_Eukarya_Species[[2]]
ggsave("Barras_Eukarya_Species.svg", plot = Barras_Eukarya_Species[[1]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
ggsave("Barras_Eukarya_Species_10.svg", plot = Barras_Eukarya_Species[[2]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Beta_Eukarya_Species<-Beta_diversity(merge_Eukaryota , 'Species' , 'Treatment', 'bray')
ggsave("Beta_Eukarya_Species.svg", plot = Beta_Eukarya_Species, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Alpha_Eukarya_Species<-Alpha_diversity(merge_Eukaryota , 'Species' , 'Treatment')
ggsave("Alpha_Eukarya_Species.svg", plot = Alpha_Eukarya_Species, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

############Bacteria############################################################


#-----------Bacteria by Phylum 
Barras_Bacteria_Phylum <- Abundance_barras(merge_Bacteria,'Phylum' , 'Treatment', 10.0)
Barras_Bacteria_Phylum[[1]] # normal
Barras_Bacteria_Phylum[[2]] # 10%
ggsave("Barras_Bacteria_Phylum.svg", plot = Barras_Bacteria_Phylum[[1]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
ggsave("Barras_Bacteria_Phylum_10.svg", plot = Barras_Bacteria_Phylum[[2]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Beta_Bacteria_Phylum<-Beta_diversity(merge_Bacteria , 'Phylum' , 'Treatment', 'bray')
ggsave("Beta_Bacteria_Phylum.svg", plot = Beta_Bacteria_Phylum, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Alpha_Bacteria_Phylum<-Alpha_diversity(merge_Bacteria , 'Phylum' , 'Treatment')
ggsave("Alpha_Bacteria_Phylum.svg", plot = Alpha_Bacteria_Phylum, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

#-----------Bacteria by Family 
Barras_Bacteria_Family <- Abundance_barras(merge_Bacteria,'Family' , 'Treatment', 10.0)
Barras_Bacteria_Family[[1]] # normal
Barras_Bacteria_Family[[2]] # 10%
ggsave("Barras_Bacteria_Family.svg", plot = Barras_Bacteria_Family[[1]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
ggsave("Barras_Bacteria_Family_10.svg", plot = Barras_Bacteria_Family[[2]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Beta_Bacteria_Family<-Beta_diversity(merge_Bacteria , 'Family' , 'Treatment', 'bray')
ggsave("Beta_Bacteria_Family.svg", plot = Beta_Bacteria_Family, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Alpha_Bacteria_Family<-Alpha_diversity(merge_Bacteria , 'Family' , 'Treatment')
ggsave("Alpha_Bacteria_Family.svg", plot = Alpha_Bacteria_Family, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

#-----------Bacteria by Genus
Barras_Bacteria_Genus <- Abundance_barras(merge_Bacteria,'Genus' , 'Treatment', 10.0)
Barras_Bacteria_Genus[[1]] # normal
Barras_Bacteria_Genus[[2]] # 10%
ggsave("Barras_Bacteria_Genus.svg", plot = Barras_Bacteria_Genus[[1]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
ggsave("Barras_Bacteria_Genus_10.svg", plot = Barras_Bacteria_Genus[[2]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Beta_Bacteria_Genus<-Beta_diversity(merge_Bacteria , 'Genus' , 'Treatment', 'bray')
ggsave("Beta_Bacteria_Genus.svg", plot = Beta_Bacteria_Genus, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Alpha_Bacteria_Genus<-Alpha_diversity(merge_Bacteria , 'Genus' , 'Treatment')
ggsave("Alpha_Bacteria_Genus.svg", plot = Alpha_Bacteria_Genus, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

#-----------Bacteria by Species 

Barras_Bacteria_Species <- Abundance_barras(merge_Bacteria,'Species' , 'Treatment', 10.0)
Barras_Bacteria_Species[[1]] 
Barras_Bacteria_Species[[2]]
ggsave("Barras_Bacteria_Species.svg", plot = Barras_Bacteria_Species[[1]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
ggsave("Barras_Bacteria_Species_10.svg", plot = Barras_Bacteria_Species[[2]], path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Beta_Bacteria_Species<-Beta_diversity(merge_Bacteria , 'Species' , 'Treatment', 'bray')
ggsave("Beta_Bacteria_Species.svg", plot = Beta_Bacteria_Species, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")

Alpha_Bacteria_Species<-Alpha_diversity(merge_Bacteria , 'Species' , 'Treatment')
ggsave("Alpha_Bacteria_Species.svg", plot = Alpha_Bacteria_Species, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
