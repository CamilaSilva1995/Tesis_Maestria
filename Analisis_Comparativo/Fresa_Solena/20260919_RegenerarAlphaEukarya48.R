## Figura 4.8 - Diversidad alfa de eucariotas por nivel taxonomico (Img/cap3/Alpha_Eukarya.png)
## Regenera la figura en alta resolucion (300 dpi) con etiquetas en espanol.
## La version anterior estaba a 960x540 px.

library("phyloseq")
library("ggplot2")
library("cowplot")

setwd("/home/camila/GIT/Tesis_Maestria/Data/fresa_solena/Data1")
outpath <- "/home/camila/GIT/Tesis_Maestria/latex/Img/cap3"

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

## Tratamiento en espanol
sam <- sample_data(fresa_kraken_fil)
sam$Treatment <- factor(sam$Treatment, levels = c("healthy", "wilted"),
                        labels = c("Saludable", "No saludable"))
sample_data(fresa_kraken_fil) <- sam

merge_Eukaryota <- subset_taxa(fresa_kraken_fil, Kingdom == "Eukaryota")

## Equivalente a Alpha_diversity() de 20230227_Funciones&Graficas.R
panel_alpha <- function(phy, tax) {
  glom <- tax_glom(phy, taxrank = tax)
  plot_richness(physeq = glom,
                measures = c("Observed", "Chao1", "Shannon", "simpson"),
                x = "Treatment", color = "Treatment") +
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

niveles <- c("Phylum", "Family", "Genus", "Species")
paneles <- lapply(niveles, function(n) panel_alpha(merge_Eukaryota, n) + guides(color = "none"))

etiquetas <- c("A) Filo", "B) Familia", "C) Género", "D) Especie")
grid_paneles <- plot_grid(plotlist = paneles, labels = etiquetas,
                          label_size = 11, hjust = 0, label_x = 0.02, ncol = 2)

## Leyenda compartida
leyenda <- get_plot_component(panel_alpha(merge_Eukaryota, "Phylum"), "guide-box-bottom")

ylab <- ggdraw() + draw_label("Medida de diversidad alfa", angle = 90, size = 12)

cuerpo <- plot_grid(ylab, grid_paneles, ncol = 2, rel_widths = c(0.035, 1))
figura <- plot_grid(cuerpo, leyenda, ncol = 1, rel_heights = c(1, 0.07)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("Alpha_Eukarya.png", plot = figura, path = outpath,
       width = 30, height = 19, dpi = 300, units = "cm", bg = "white")

cat("Listo 4.8\n")

ggsave("Alpha_Eukarya.png", plot = figura, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img",
       width = 30, height = 19, dpi = 300, units = "cm", bg = "white")
