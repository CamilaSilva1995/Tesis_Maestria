## Figura 4.1 - Diversidad alfa antes y despues del filtro de calidad
## (Img/cap3/AlphaDiversity_CrudosFiltrados.png)
## Regenera la figura en alta resolucion (300 dpi) con etiquetas en espanol.
## La version anterior estaba a 960x349 px. El script original no quedo en el
## repositorio; esta version reproduce el analisis de 01_Exploracion.Rmd
## (lineas 268 y 275), que grafica plot_richness sobre el objeto crudo y sobre
## el filtrado por calidad.

library("phyloseq")
library("ggplot2")
library("cowplot")

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
outpath_analisis <- "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img"
outpath_tesis <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap3"

fresa_kraken <- import_biom("fresa_kraken.biom")
colnames(fresa_kraken@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data, 4, 100)
colnames(fresa_kraken@otu_table@.Data) <- substr(colnames(fresa_kraken@otu_table@.Data), 1, 6)
metadata_fresa <- read.csv2("metadata.csv", header = FALSE, row.names = 1, sep = ",")
fresa_kraken@sam_data <- sample_data(metadata_fresa)
fresa_kraken@sam_data$Sample <- row.names(fresa_kraken@sam_data)
colnames(fresa_kraken@sam_data) <- c('Treatment', 'Samples')

## Tratamiento en espanol
sam <- sample_data(fresa_kraken)
sam$Treatment <- factor(sam$Treatment, levels = c("healthy", "wilted"),
                        labels = c("Saludable", "No saludable"))
sample_data(fresa_kraken) <- sam

samples_to_remove <- c("MP2079", "MP2080", "MP2088", "MP2109", "MP2137")
fresa_kraken_fil <- prune_samples(!(sample_names(fresa_kraken) %in% samples_to_remove), fresa_kraken)

cat("muestras crudas:", nsamples(fresa_kraken),
    " filtradas:", nsamples(fresa_kraken_fil), "\n")

panel_alpha <- function(phy) {
  p <- plot_richness(physeq = phy,
                     measures = c("Observed", "Chao1", "Shannon", "simpson"),
                     x = "Treatment", color = "Treatment")
  ## plot_richness rotula la medida "Observed" en ingles
  p$data$variable <- factor(p$data$variable,
                            levels = c("Observed", "Chao1", "Shannon", "Simpson"),
                            labels = c("Observado", "Chao1", "Shannon", "Simpson"))
  p +
    labs(x = NULL, y = NULL) +
    scale_color_discrete(name = "Tratamiento") +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 10),
          text = element_text(size = 11),
          strip.text = element_text(size = 9),
          axis.text.x = element_text(angle = 45, size = 8, hjust = 1),
          axis.text.y = element_text(size = 7),
          plot.margin = margin(20, 6, 42, 6))
}

panelA <- panel_alpha(fresa_kraken)     # datos crudos
panelB <- panel_alpha(fresa_kraken_fil) # tras el filtro de calidad

paneles <- plot_grid(panelA + guides(color = "none"),
                     panelB + guides(color = "none"),
                     labels = c("A) Datos crudos", "B) Tras el filtro de calidad"),
                     label_size = 11, hjust = 0, label_x = 0.01, ncol = 2)

leyenda <- get_plot_component(panelA, "guide-box-bottom")
ylab <- ggdraw() + draw_label("Medida de diversidad alfa", angle = 90, size = 12)

cuerpo <- plot_grid(ylab, paneles, ncol = 2, rel_widths = c(0.035, 1))
figura <- plot_grid(cuerpo, leyenda, ncol = 1, rel_heights = c(1, 0.08)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("AlphaDiversity_CrudosFiltrados.png", plot = figura, path = outpath_analisis,
       width = 30, height = 15, dpi = 300, units = "cm", bg = "white")
ggsave("AlphaDiversity_CrudosFiltrados.png", plot = figura, path = outpath_tesis,
       width = 30, height = 15, dpi = 300, units = "cm", bg = "white")

cat("Listo 4.1\n")
