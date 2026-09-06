## Regenera en alta resolucion (300 dpi) la figura de barras de abundancia por Filo
## (Figura 4.5 de la tesis, Img/cap3/BarrasFilo_blackandwhite.png), con etiquetas en
## espanol. El script original que genero esta figura no quedo en el repositorio
## (se genero de forma interactiva); esta version usa la misma funcion
## Abundance_barras() de 20230227_Funciones&Graficas.R y la misma paleta de color
## por Filo (colorRampPalette(brewer.pal(8,"Dark2"))) usada en el resto del proyecto
## (ver 02_ExploracionSubconjuntos.Rmd, 07_RedesCoocurrencia.Rmd).

library("phyloseq")
library("ggplot2")
library("RColorBrewer")

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
outpath_analisis <- "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"
outpath_tesis <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap3"

### Cargado de datos originales (igual que 20230227_Funciones&Graficas.R)
fresa_kraken <- import_biom("fresa_kraken.biom")
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data, 4, 100)
colnames(fresa_kraken@otu_table@.Data) <- substr(colnames(fresa_kraken@otu_table@.Data), 1, 6)
metadata_fresa <- read.csv2("metadata.csv", header = FALSE, row.names = 1, sep = ",")
fresa_kraken@sam_data <- sample_data(metadata_fresa)
fresa_kraken@sam_data$Sample <- row.names(fresa_kraken@sam_data)
colnames(fresa_kraken@sam_data) <- c('Treatment', 'Samples')

samples_to_remove <- c("MP2079", "MP2080", "MP2088", "MP2109", "MP2137")
fresa_kraken_fil <- prune_samples(!(sample_names(fresa_kraken) %in% samples_to_remove), fresa_kraken)

glomToGraph <- function(phy, tax) {
  glom <- tax_glom(phy, taxrank = tax)
  percentages <- transform_sample_counts(glom, function(x) x * 100 / sum(x))
  percentages_df <- psmelt(percentages)
  return(list(glom, percentages, percentages_df))
}

## Version en espanol de Abundance_barras(), con la paleta de color estandar del
## proyecto para el relleno por Filo y borde blanco/negro por Tratamiento.
Abundance_barras_es <- function(phy, tax, attribute, abundance_percentage) {
  Data <- glomToGraph(phy, tax)
  percentages_df <- Data[[3]]
  percentages_df[[tax]] <- as.factor(percentages_df[[tax]])
  tax_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(levels(percentages_df[[tax]])))
  plot_barras <- ggplot(data = percentages_df, aes_string(x = 'Sample', y = 'Abundance', fill = tax, color = attribute)) +
    scale_fill_manual(values = tax_colors) +
    scale_colour_manual(name = "Tratamiento", values = c('white', 'black'), labels = c("Saludable", "No saludable")) +
    geom_bar(aes(), stat = "identity", position = "stack") +
    labs(title = "Abundancia", x = 'Muestra', y = 'Abundancia', fill = 'Filo') +
    theme(legend.key.size = unit(0.2, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 6),
          text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, size = 5, hjust = 1, vjust = 0.5))
  return(plot_barras)
}

filo_plot_blackandwhite_es <- Abundance_barras_es(fresa_kraken_fil, 'Phylum', 'Treatment', 10.0)

ggsave("BarrasFilo_blackandwhite.png", plot = filo_plot_blackandwhite_es, path = outpath_analisis,
       width = 30, height = 15, dpi = 300, units = "cm")
ggsave("BarrasFilo_blackandwhite.png", plot = filo_plot_blackandwhite_es, path = outpath_tesis,
       width = 30, height = 15, dpi = 300, units = "cm")
