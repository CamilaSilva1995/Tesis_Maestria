## Figura 4.6 - Abundancia relativa de eucariotas por nivel taxonomico (Img/cap3/Barras_Eukarya10.png)
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

## Tratamiento en espanol (healthy -> Saludable, wilted -> No saludable)
sam <- sample_data(fresa_kraken_fil)
sam$Treatment <- factor(sam$Treatment, levels = c("healthy", "wilted"),
                        labels = c("Saludable", "No saludable"))
sample_data(fresa_kraken_fil) <- sam

merge_Eukaryota <- subset_taxa(fresa_kraken_fil, Kingdom == "Eukaryota")

## Equivalente al plot_percentages (variante [[2]]) de Abundance_barras():
## agrupa en "Otros" los taxones por debajo del umbral de abundancia.
panel_barras <- function(phy, tax, umbral) {
  glom <- tax_glom(phy, taxrank = tax)
  percentages <- transform_sample_counts(glom, function(x) x * 100 / sum(x))
  df <- psmelt(percentages)
  df$tax <- df[[tax]]
  df$tax[df$Abundance < umbral] <- "Otros"
  df$tax <- as.factor(df$tax)
  ggplot(data = df, aes(x = Sample, y = Abundance, fill = tax, color = Treatment)) +
    scale_colour_manual(name = "Tratamiento", values = c('white', 'black')) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = NULL, y = NULL, fill = NULL) +
    theme(legend.key.size = unit(0.22, "cm"),
          legend.key.width = unit(0.22, "cm"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 6.5),
          text = element_text(size = 11),
          axis.text.x = element_text(angle = 90, size = 4, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 8),
          plot.margin = margin(20, 6, 4, 6))
}

niveles <- list(c("Phylum", "Filo"), c("Family", "Familia"),
                c("Genus", "Genero"), c("Species", "Especie"))

paneles <- lapply(niveles, function(n) {
  p <- panel_barras(merge_Eukaryota, n[1], 10.0)
  ## cada panel conserva su leyenda de taxones; la de Tratamiento es compartida
  p + guides(colour = "none")
})

etiquetas <- c("A) Filo", "B) Familia", "C) Género", "D) Especie")
grid_paneles <- plot_grid(plotlist = paneles, labels = etiquetas,
                          label_size = 11, hjust = 0, label_x = 0.02, ncol = 2)

## Leyenda compartida de Tratamiento
p_trat <- panel_barras(merge_Eukaryota, "Phylum", 10.0) + guides(fill = "none")
leyenda_trat <- get_plot_component(p_trat, "guide-box-bottom")

## Etiquetas de ejes compartidas
ylab <- ggdraw() + draw_label("Abundancia", angle = 90, size = 12)
xlab <- ggdraw() + draw_label("Muestra", size = 12)

cuerpo <- plot_grid(ylab, grid_paneles, ncol = 2, rel_widths = c(0.03, 1))
figura <- plot_grid(cuerpo, xlab, leyenda_trat, ncol = 1,
                    rel_heights = c(1, 0.04, 0.06)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("Barras_Eukarya10.png", plot = figura, path = outpath,
       width = 30, height = 17, dpi = 300, units = "cm", bg = "white")

cat("Listo 4.6\n")

ggsave("Barras_Eukarya10.png", plot = figura, path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img",
       width = 30, height = 17, dpi = 300, units = "cm", bg = "white")
