### Crear archivos .biom 
# en la carpeta solena, delosmeta datos tomamos solo los de tomate
# grep Tomate fastp_metadat.csv -> vemos solo los de tomate
# grep Tomate fastp_metadat.csv | cut -d "," -f 1 > tomate.txt -> se filtra solo por tomate, la primera fila "nombre" y lo guardamos en un .txt
# luego entrando con vi modificamos el archivo para que quede como script de bash y se pueda correr el kraken-biom 
# (para reemplazar caracteres en vi se usa para reemplazar un caracter por estacio, por ejemplo: ":s/./ /" o los  saltos de linea por estacios ":s/\n/ /" y si lo que se quiere reemplazar son varias cosas en la misma linea se termina el comando con una g,como por ejemplo: ":s/./ /g")
# kraken-biom vcXtM2TOMPOOB050180.report iYASM2TOMPOOA050195.report XR8LM2TOMPOOB051125.report TECcM2TOMPOOA051310.report FvY2M2TOMPOOB051329.report abKIM2TOMPOOA051320.report cuu6M2TOMPOOB051315.report s9OmM2TOMPOOB051285.report Nm8MM2TOMPOOA051300.report IGpcM2TOMPOOB051290.report a50wM2TOMPOOB051350.report OeS0M2TOMPOOA051355.report dSz8M2TOMPOOA051340.report aFJgM2TOMPOOM090325.report cSObM2TOMPOOM080315.report wQqUM2TOMPOOA080335.report 0n7bM2TOMPOOA100365.report 1qGgM2TOMPOOM080340.report 9JUZM2TOMPOOB100355.report AYtCM2TOMPOOA100360.report DYZgM2TOMPOOB100370.report nG3jM2TOMPOOA051345.report --fmt json -o tomatecamila.biom
# para correr el kraken-biom es necesario entrar en ambiente de conda metagenomics, con: conda activate metagenomic
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
merged_metagenomes <- import_biom("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/tomatecamila.biom")
class(merged_metagenomes)

### metagenoma rizosfera del tomate

# vemos cuantos reads tenemos de cada muestra
sample_sums(merged_metagenomes)

## vemos la tabla de clasificacion taxonomica
View(merged_metagenomes@tax_table@.Data)
# se remueven catarteres inecesarios
merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
# renombramos las columnas 
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# para saber cuantos phylum diferentes hay
unique(merged_metagenomes@tax_table@.Data[,"Phylum"])


## vemos la tabla de abundancia taxonomica 
View(merged_metagenomes@otu_table@.Data)

summary(merged_metagenomes@otu_table@.Data)

# visualizar la diversidad (Alpha)
pdf("alphadiversity_tomate_plot.pdf") 
p <- plot_richness(physeq = merged_metagenomes, 
                   measures = c("Observed","Chao1","Shannon")) 
# cierra el archivo pdf
dev.off() 
ggsave("alphadiversity_tomate_plot.pdf", plot = p)

# visualizar la diversidad (Beta)
# trasformamos latabla deabuncandia a porcentajes
percentages <- transform_sample_counts(merged_metagenomes, function(x) x*100 / sum(x) )

meta_ord <- ordinate(physeq = percentages, method = "NMDS", 
                     distance = "bray")

pdf("betadiversity_tomate_plot.pdf") 
p2 <- plot_ordination(physeq = percentages, ordination = meta_ord) 
# cierra el archivo pdf
dev.off() 
ggsave("betadiversity_tomate_plot.pdf", plot = p2)




###################################################################
### Repetir todo con el maiz y el chile ###
### Ravisar el mismo procedimiento con el paquete "VEGAN" 


