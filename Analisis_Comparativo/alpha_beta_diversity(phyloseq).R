### Crear archivos .biom 
## https://github.com/carpentries-incubator/metagenomics
## https://carpentries-incubator.github.io/metagenomics/
# en la carpeta solena, de los meta datos tomamos solo los de tomate
# grep Tomate fastp_metadat.csv -> vemos solo los de tomate
# grep Tomate fastp_metadat.csv | cut -d "," -f 1 > tomate.txt -> se filtra solo por tomate, la primera fila "nombre" y lo guardamos en un .txt
# luego entrando con vi modificamos el archivo para que quede como script de bash y se pueda correr el kraken-biom 
# /home/betterlab/GIT/Tesis_Maestria/Analisis_Comparativo
# (para reemplazar caracteres en vi se usa para reemplazar un caracter por estacio, por ejemplo: ":s/./ /" o los  saltos de linea por estacios ":s/\n/ /" y si lo que se quiere reemplazar son varias cosas en la misma linea se termina el comando con una g,como por ejemplo: ":s/./ /g")
# kraken-biom vcXtM2TOMPOOB050180.report iYASM2TOMPOOA050195.report XR8LM2TOMPOOB051125.report TECcM2TOMPOOA051310.report FvY2M2TOMPOOB051329.report abKIM2TOMPOOA051320.report cuu6M2TOMPOOB051315.report s9OmM2TOMPOOB051285.report Nm8MM2TOMPOOA051300.report IGpcM2TOMPOOB051290.report a50wM2TOMPOOB051350.report OeS0M2TOMPOOA051355.report dSz8M2TOMPOOA051340.report aFJgM2TOMPOOM090325.report cSObM2TOMPOOM080315.report wQqUM2TOMPOOA080335.report 0n7bM2TOMPOOA100365.report 1qGgM2TOMPOOM080340.report 9JUZM2TOMPOOB100355.report AYtCM2TOMPOOA100360.report DYZgM2TOMPOOB100370.report nG3jM2TOMPOOA051345.report --fmt json -o tomatecamila.biom
# para correr el kraken-biom es necesario entrar en ambiente de conda metagenomics, con: conda activate metagenomics
# se corre "sh tomate.sh" y asi se crea un archivo.biom "tomatecamila.biom"
# luego de esto debemos ver si quedo bien en R

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  +     install.packages("BiocManager")
#BiocManager::install("phyloseq") 

#install.packages(c("ggplot2", "readr", "patchwork"))

library("phyloseq")
library("ggplot2")
library("readr")
library("patchwork")
getwd()
setwd("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo")
merged_metagenomes_tomate <- import_biom("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/tomatecamila.biom")
class(merged_metagenomes_tomate)

### metagenoma rizosfera del tomate

# vemos cuantos reads tenemos de cada muestra
sample_sums(merged_metagenomes_tomate)

## vemos la tabla de clasificacion taxonomica
View(merged_metagenomes_tomate@tax_table@.Data)
# se remueven catarteres inecesarios
merged_metagenomes_tomate@tax_table@.Data <- substring(merged_metagenomes_tomate@tax_table@.Data, 4)
# renombramos las columnas 
colnames(merged_metagenomes_tomate@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# para saber cuantos phylum diferentes hay
unique(merged_metagenomes_tomate@tax_table@.Data[,"Phylum"])


## vemos la tabla de abundancia taxonomica 
View(merged_metagenomes_tomate@otu_table@.Data)

summary(merged_metagenomes_tomate@otu_table@.Data)

# visualizar la diversidad (Alpha)
pdf("alphadiversity_tomate_plot.pdf") 
p_tomate <- plot_richness(physeq = merged_metagenomes_tomate, 
                          measures = c("Observed","Chao1","Shannon")) 
# cierra el archivo pdf
dev.off() 
ggsave("alphadiversity_tomate_plot.pdf", plot = p_tomate)

# visualizar la diversidad (Beta)
# trasformamos latabla deabuncandia a porcentajes
percentages_tomate <- transform_sample_counts(merged_metagenomes_tomate, function(x) x*100 / sum(x) )

meta_ord_tomate <- ordinate(physeq = percentages_tomate, method = "NMDS", 
                            distance = "bray")

pdf("betadiversity_tomate_plot.pdf") 
p2_tomate <- plot_ordination(physeq = percentages_tomate, ordination = meta_ord_tomate) 
# cierra el archivo pdf
dev.off() 
ggsave("betadiversity_tomate_plot.pdf", plot = p2_tomate)




###################################################################

merged_metagenomes_maiz <- import_biom("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/maizcamila.biom")
class(merged_metagenomes_maiz)

### metagenoma rizosfera del maiz

# vemos cuantos reads tenemos de cada muestra
sample_sums(merged_metagenomes_maiz)

## vemos la tabla de clasificacion taxonomica
View(merged_metagenomes_maiz@tax_table@.Data)
# se remueven catarteres inecesarios
merged_metagenomes_maiz@tax_table@.Data <- substring(merged_metagenomes_maiz@tax_table@.Data, 4)
# renombramos las columnas 
colnames(merged_metagenomes_maiz@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# para saber cuantos phylum diferentes hay
unique(merged_metagenomes_maiz@tax_table@.Data[,"Phylum"])


## vemos la tabla de abundancia taxonomica 
View(merged_metagenomes_maiz@otu_table@.Data)

summary(merged_metagenomes_maiz@otu_table@.Data)

# visualizar la diversidad (Alpha)
pdf("alphadiversity_maiz_plot.pdf") 
p_maiz<- plot_richness(physeq = merged_metagenomes_maiz, 
                       measures = c("Observed","Chao1","Shannon")) 
# cierra el archivo pdf
dev.off() 
ggsave("alphadiversity_maiz_plot.pdf", plot = p_maiz)

# visualizar la diversidad (Beta)
# trasformamos latabla deabuncandia a porcentajes
percentages_maiz <- transform_sample_counts(merged_metagenomes_maiz, function(x) x*100 / sum(x) )

meta_ord_maiz <- ordinate(physeq = percentages_maiz, method = "NMDS", 
                          distance = "bray")

pdf("betadiversity_maiz_plot.pdf") 
p2_maiz <- plot_ordination(physeq = percentages_maiz, ordination = meta_ord_maiz) 
# cierra el archivo pdf
dev.off() 
ggsave("betadiversity_maiz_plot.pdf", plot = p2_maiz)


###################################################################

merged_metagenomes_chile <- import_biom("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/chile.biom")
class(merged_metagenomes_chile)

### metagenoma rizosfera del chile

# vemos cuantos reads tenemos de cada muestra
sample_sums(merged_metagenomes_chile)

## vemos la tabla de clasificacion taxonomica
View(merged_metagenomes_chile@tax_table@.Data)
# se remueven catarteres inecesarios
merged_metagenomes_chile@tax_table@.Data <- substring(merged_metagenomes_chile@tax_table@.Data, 4)
# renombramos las columnas 
colnames(merged_metagenomes_chile@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# para saber cuantos phylum diferentes hay
unique(merged_metagenomes_chile@tax_table@.Data[,"Phylum"])


## vemos la tabla de abundancia taxonomica 
View(merged_metagenomes_chile@otu_table@.Data)

summary(merged_metagenomes_chile@otu_table@.Data)

# visualizar la diversidad (Alpha)
pdf("alphadiversity_chile_plot.pdf") 
p_chile <- plot_richness(physeq = merged_metagenomes_chile, 
                         measures = c("Observed","Chao1","Shannon")) 
# cierra el archivo pdf
dev.off() 
ggsave("alphadiversity_chile_plot.pdf", plot = p_chile)

# visualizar la diversidad (Beta)
# trasformamos latabla deabuncandia a porcentajes
percentages_chile <- transform_sample_counts(merged_metagenomes_chile, function(x) x*100 / sum(x) )

meta_ord_chile <- ordinate(physeq = percentages_chile, method = "NMDS", 
                           distance = "bray")

pdf("betadiversity_chile_plot.pdf") 
p2_chile <- plot_ordination(physeq = percentages_chile, ordination = meta_ord_chile) 
# cierra el archivo pdf
dev.off() 
ggsave("betadiversity_chile_plot.pdf", plot = p2_chile)

### Ravisar el mismo procedimiento con el paquete "VEGAN" 



