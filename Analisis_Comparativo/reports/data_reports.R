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

# para explorara que Phylum tenemos 
unique(merged_metagenomes@tax_table@.Data[,"Phylum"])

#open a pdf to save the plot
pdf("richness_plot.pdf") 
plot_richness(physeq = merged_metagenomes, 
                measures = c("Observed","Chao1","Shannon")) 
# Close the pdf file
dev.off() 



# metadatos... tres tablas... maiz/chile/tomate





