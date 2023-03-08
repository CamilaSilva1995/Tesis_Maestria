## CREAR UNA FUNCION QUE GRAFIQUE TODAS LAS OPCIONES

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
#-----------------------------------------------------------------
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
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold"),
          text = element_text(size=12),
          axis.text.x = element_text(angle=90, size=12, hjust=1, vjust=0.5))
  
  percentages_df$tax<-percentages_df[,ncol(percentages_df)]
  percentages_df$tax[percentages_df$Abundance < abundance_percentage] <- "abundance_percentage"
  percentages_df$tax <- as.factor(percentages_df$tax)
  plot_percentages <- ggplot(data=percentages_df, aes_string(x='Sample', y='Abundance', fill='tax' ,color=attribute))+
    scale_colour_manual(values=c('white','black')) +
    geom_bar(aes(), stat="identity", position="stack") +
    labs(title = "Abundance", x='Sample', y='Abundance', color = tax) +
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold"),
          text = element_text(size=12),
          axis.text.x = element_text(angle=90, size=12, hjust=1, vjust=0.5))
  return(list(plot_barras,plot_percentages))
}

Barras <- Abundance_barras(merge_Bacteria, 'Genus' , 'Treatment', 10.0)
Barras[1] # normal
Barras[2]





#-----------------------------------------
#GRaficar betadiversity
Beta_diversity <- function(phy,tax,attribute,distance){
  ##Podemos llamar a la funcion que crea los glom ??
  Data <- glomToGraph(phy,tax)
  glom <- Data[[1]]
  percentages <- Data[[2]]
  percentages_df <- Data[[3]]
  ## Beta diversidad 
  meta_ord <- ordinate(physeq = percentages, method = "NMDS", distance = distance) 
  plot_beta <- plot_ordination(physeq = percentages, ordination = meta_ord, color = attribute) +
    geom_text(mapping = aes(label = colnames(phy@otu_table@.Data)), size = 3, vjust = 1.5)
  return(plot_beta)
}

Alpha_diversity <- function(phy,tax,attribute){
  ##Podemos llamar a la funcion que crea los glom ??
  Data <- glomToGraph(phy,tax)
  glom <- Data[[1]]
  percentages <- Data[[2]]
  percentages_df <- Data[[3]]
  ## Alfa diversidad
  plot_alpha <- plot_richness(physeq = glom, measures = c("Observed","Chao1","Shannon","simpson"),x = attribute, color = attribute) 
  #ggsave("DiversidadAlfa_tax_fil.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
  return(plot_alpha)
}

###################################################################3
#### main program 
####################################################################3

#-----------Eukarya by Phylum 

#glomToGraph(merge_Eukaryota,'Phylum') 
Barras_Phylum <- Abundance_barras(merge_Eukaryota,'Phylum' , 'Treatment', 10.0)
Barras_Phylum[1] # normal
Barras_Phylum[2]
# ggsave("Barras_Phylum_Eukaryotanew.png", plot = last_plot(), path = outpath , width = 30, height = 15, dpi = 300, units = "cm")
Beta_diversity(merge_Eukaryota , 'Phylum' , 'Treatment', 'bray')

Alpha_diversity(merge_Eukaryota , 'Phylum' , 'Treatment')


#-----------Eukarya by Species 

Barras_Phylum <- Abundance_barras(merge_Eukaryota,'Species' , 'Treatment', 10.0)
Barras_Phylum[1] 
Barras_Phylum[2]

#-----------Bacteria by Species

Barras_Phylum <- Abundance_barras(merge_Bacteria,'Species' , 'Treatment', 10.0)
Barras_Phylum[1] 
Barras_Phylum[2]








