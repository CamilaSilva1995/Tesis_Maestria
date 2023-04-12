## GRAFICAS DE PROFUNDIDAD Y NORMALIZACION 

library("phyloseq")
library("ggplot2")
library("vegan")
#library("BiodiversityR") 
library("RColorBrewer")
library("stringi")
library("dplyr")
library("plyr")

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
outpath = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"
### Cargado de datos 
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


## Queremos ver la profundidad de las muestras
dprof <- data.frame(Samples = colnames(fresa_kraken_fil@otu_table@.Data),
                    Reads = sample_sums(fresa_kraken_fil),
                    Treatment = fresa_kraken_fil@sam_data@.Data[[1]])

mu_Samples <- mean(dprof$Reads) #ddply(dprof, "Treatment", summarise, grp.mean=mean(Reads))#

ggplot(data = dprof, mapping = aes(x = Samples, y = Reads))+
  geom_bar(stat = "identity", aes( fill = Treatment)) +
  geom_hline(yintercept=mu_Samples, color = "black", linetype="dashed", size=1.5) +
  scale_fill_manual(values = c("cyan3","darkmagenta")) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# La línea marca la media de la profundidad de las muestras.
# 
# Esto deja en claro que los datos necesitan normalización. McMurdi et al. 
# encontró una gran metodología para normalizar los datos. En su artículo Waste Not, 
# Want Not: Why Rarefying Microbiome Data Is Inadmisible, hacen una interesante
# discusión sobre este tema constante en este tipo de análisis. Usaremos su metodología:
#

edgeRnorm = function(phy, ...){
  require("edgeR")
  require("phyloseq")
  if (!taxa_are_rows(phy)) {
    phy <- t(phy)
  }
  x = as(otu_table(phy), "matrix")
  x = x + 1
  y = edgeR::DGEList(counts = x, remove.zeros = TRUE)
  z = edgeR::calcNormFactors(y, ...)
  if (!all(is.finite(z$samples$norm.factors))) {
    stop("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors")
  }
  return(z)
}

z<- edgeRnorm(cbact, method = "TMM")

nor.cb <- merge_phyloseq(otu_table(z@.Data[[1]], taxa_are_rows = TRUE),
                         tax_table(cbact@tax_table@.Data),
                         cbact@sam_data)
rm(z)





























































