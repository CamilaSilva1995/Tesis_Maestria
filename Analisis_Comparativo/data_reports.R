#https://carpentries-incubator.github.io/metagenomics/07-Diversity-tackled-with-R/index.html
getwd()
setwd("/home/betterlab/GIT/Tesis_Maestria/Analisis_Comparativo/reports")


if (!requireNamespace("BiocManager", quietly = TRUE))
  +     install.packages("BiocManager")
BiocManager::install("phyloseq") 

install.packages(c("ggplot2", "readr", "patchwork"))

library("phyloseq")
library("ggplot2")
library("readr")
library("patchwork")
getwd()
merged_metagenomes <- import_biom("data_reports.biom")
class(merged_metagenomes)


View(merged_metagenomes@tax_table@.Data)
merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
View(merged_metagenomes@tax_table@.Data)
# ponemos los nombres a lascolumnas 
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
View(merged_metagenomes@tax_table@.Data)

# para explorar a que Phylum tenemos 
unique(merged_metagenomes@tax_table@.Data[,"Phylum"])
# cuantos OTUS estan clasificados en filo Firmicutes
sum(merged_metagenomes@tax_table@.Data[,"Phylum"] == "Firmicutes")
# cuantos son "Proteobacteria"
sum(merged_metagenomes@tax_table@.Data[,"Phylum"] == "Proteobacteria")




merged_metagenomes <- subset_taxa(merged_metagenomes, Kingdom == "Bacteria")
merged_metagenomes

sample_sums(merged_metagenomes)
# tenemos una muestra de 1048070 reads
summary(merged_metagenomes@otu_table@.Data)


# queremos visualizar la diversidad (Alpha)
pdf("richness_plot.pdf") 
p <- plot_richness(physeq = merged_metagenomes, 
                measures = c("Observed","Chao1","Shannon")) 
# Close the pdf file
dev.off() 
ggsave("richness_plot.pdf", plot = p)


#######################################################################
# metadatos... tres tablas... maiz/chile/tomate
# 
library("phyloseq")
library("ggplot2")
library("edgeR")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("stringr")
library("tidyverse")
library("vegan")

solena <- import_biom("/home/camila/GIT/Tesis_Maestria/solena/biom/solena_maiz.biom")
solena

solena@tax_table@.Data <- substring(solena@tax_table@.Data, 4)
colnames(solena@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(solena@otu_table@.Data) <- substring(sample_names(solena), 3)

colnames(solena@otu_table@.Data) <- substring(sample_names(solena), 15)

# Cargar los metadatos de solena en una tabla.
met.sol <- read.csv2("/home/camila/GIT/Tesis_Maestria/solena/fastp_metadat.csv",
                     header =  TRUE, row.names = 1, sep = ",")
rownames(met.sol) <- substring(rownames(met.sol),first = 13)
# se trabaja solo con una muestra de los datos
met.sol <- sample_data(met.sol)

solena<- merge_phyloseq(solena, met.sol)
solena

# Profundidad de los datos
dep.sol <- data.frame(Samples = as.factor(colnames(solena@otu_table@.Data)),
                      Reads = sample_sums(solena),
                      Cultivo = as.factor(solena@sam_data@.Data[[3]]))
