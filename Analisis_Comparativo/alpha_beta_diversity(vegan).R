library("phyloseq")
library("ggplot2")
library("readr")
library("patchwork")
library("vegan")

getwd()
setwd("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo")
merged_metagenomes_tomate <- import_biom("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/tomatecamila.biom")
class(merged_metagenomes_tomate)

abundance_table <- merged_metagenomes_tomate@otu_table@.Data
abundance_table <- t(abundance_table)
H<-diversity(abundance_table, "shannon")
df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon", length(H)))
