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



# metadatos... tres tablas... maiz/chile/tomate
# 




