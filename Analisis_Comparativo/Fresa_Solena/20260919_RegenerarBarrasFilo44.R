## Regenera en alta resolucion (300 dpi) la figura de barras de abundancia por Filo
## con dos paneles (Figura 4.4 de la tesis, Img/cap3/BarrasFilo.png):
##   A) abundancias absolutas  B) abundancias relativas
## con etiquetas en espanol y leyenda compartida.
##
## La version anterior estaba a 960x540 px, con las etiquetas de muestra
## ilegibles. El script original no quedo en el repositorio (se genero de forma
## interactiva); esta version reproduce el mismo analisis a partir del codigo de
## 02_ExploracionSubconjuntos.Rmd (chunks absolute_plot y relative_plot) y usa la
## misma paleta por Filo, colorRampPalette(brewer.pal(8,"Dark2")), empleada en el
## resto del proyecto y en la Figura 4.5.

library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("cowplot")

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
outpath_analisis <- "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"
outpath_tesis <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap3"

### Cargado de datos originales (igual que 02_ExploracionSubconjuntos.Rmd)
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
percentages_fil <- transform_sample_counts(fresa_kraken_fil, function(x) x * 100 / sum(x))

## Panel A: abundancias absolutas (conteo de lecturas por filo)
absolute_glom_phylum <- tax_glom(fresa_kraken_fil, taxrank = "Phylum")
absolute_df_phylum <- psmelt(absolute_glom_phylum)
absolute_df_phylum$Phylum <- as.factor(absolute_df_phylum$Phylum)

## Panel B: abundancias relativas (porcentaje por muestra)
percentages_glom_phylum <- tax_glom(percentages_fil, taxrank = "Phylum")
percentages_df_phylum <- psmelt(percentages_glom_phylum)
percentages_df_phylum$Phylum <- as.factor(percentages_df_phylum$Phylum)

phylum_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(levels(absolute_df_phylum$Phylum)))

barras_theme <- theme(
  legend.key.size = unit(0.35, "cm"),
  legend.key.width = unit(0.35, "cm"),
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.title = element_text(size = 11, face = "bold"),
  legend.text = element_text(size = 9),
  text = element_text(size = 12),
  axis.text.x = element_text(angle = 90, size = 6, hjust = 1, vjust = 0.5),
  axis.text.y = element_text(size = 9)
)

hacer_barras <- function(df) {
  ggplot(data = df, aes(x = Sample, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = phylum_colors) +
    labs(x = "Muestra", y = "Abundancia", fill = "Filo") +
    guides(fill = guide_legend(nrow = 11)) +
    barras_theme
}

absolute_plot <- hacer_barras(absolute_df_phylum)
relative_plot <- hacer_barras(percentages_df_phylum)

## Leyenda compartida. Se usa get_plot_component en lugar de get_legend porque
## en cowplot >= 1.1.3 este ultimo devuelve un guide-box vacio.
leyenda <- get_plot_component(absolute_plot, "guide-box-bottom")

paneles <- plot_grid(
  absolute_plot + theme(legend.position = "none",
                        plot.margin = margin(5, 10, 5, 5)),
  relative_plot + theme(legend.position = "none",
                        plot.margin = margin(5, 5, 5, 10)),
  labels = c('A)', 'B)'), label_size = 14, ncol = 2
)

figura <- plot_grid(paneles, leyenda, ncol = 1, rel_heights = c(1, 0.5)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("BarrasFilo.png", plot = figura, path = outpath_analisis,
       width = 30, height = 19, dpi = 300, units = "cm", bg = "white")
ggsave("BarrasFilo.png", plot = figura, path = outpath_tesis,
       width = 30, height = 19, dpi = 300, units = "cm", bg = "white")
