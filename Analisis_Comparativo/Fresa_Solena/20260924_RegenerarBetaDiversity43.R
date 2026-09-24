## Figura 4.3 - Diversidad beta (NMDS) con cinco metricas de distancia
## (Img/cap3/BetaDiversity.png)
## Regenera la figura en alta resolucion (300 dpi) con etiquetas en espanol.
##
## La version anterior estaba a 960x520 px. El script original no quedo en el
## repositorio; esta version reproduce las ordenaciones de 01_Exploracion.Rmd
## (lineas 372-426). Se omiten las etiquetas de nombre de muestra que llevaba
## la version anterior (geom_text sobre cada punto): con 53 muestras se
## encimaban unas sobre otras, eran ilegibles y tapaban los propios puntos que
## la figura busca comparar.

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

samples_to_remove <- c("MP2079", "MP2080", "MP2088", "MP2109", "MP2137")
fresa_kraken_fil <- prune_samples(!(sample_names(fresa_kraken) %in% samples_to_remove), fresa_kraken)

sam <- sample_data(fresa_kraken_fil)
sam$Treatment <- factor(sam$Treatment, levels = c("healthy", "wilted"),
                        labels = c("Saludable", "No saludable"))
sample_data(fresa_kraken_fil) <- sam

percentages_fil <- transform_sample_counts(fresa_kraken_fil, function(x) x * 100 / sum(x))

## NMDS es estocastico: la semilla fija hace reproducible la figura
set.seed(20260924)

## Las ordenaciones NMDS tardan varios minutos; se cachean en disco para poder
## reajustar la figura sin recalcularlas. Borrar el .rds fuerza el recalculo.
cache_ord <- file.path(outpath_analisis, "BetaDiversity_ordinaciones.rds")

panel_beta <- function(distancia) {
  ord <- ordenaciones[[distancia]]
  cat("  ", distancia, "- stress:", round(ord$stress, 4), "\n")
  p <- plot_ordination(physeq = percentages_fil, ordination = ord, color = "Treatment") +
    geom_point(size = 2.2) +
    scale_color_discrete(name = "Tratamiento") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 10),
          text = element_text(size = 11),
          axis.text = element_text(size = 8),
          plot.margin = margin(18, 8, 6, 6))
  list(plot = p, stress = ord$stress)
}

metricas <- c("bray", "jaccard", "euclidean", "manhattan", "jsd")
nombres  <- c("Bray-Curtis", "Jaccard", "Euclidiana", "Manhattan", "Jensen-Shannon")

if (file.exists(cache_ord)) {
  cat("Reutilizando ordenaciones cacheadas\n")
  ordenaciones <- readRDS(cache_ord)
} else {
  cat("Calculando ordenaciones NMDS (esto tarda varios minutos):\n")
  ordenaciones <- lapply(metricas, function(d)
    ordinate(physeq = percentages_fil, method = "NMDS", distance = d))
  names(ordenaciones) <- metricas
  saveRDS(ordenaciones, cache_ord)
}

resultados <- lapply(metricas, panel_beta)

etiquetas <- paste0(LETTERS[1:5], ") ", nombres,
                    " (estres = ", sprintf("%.3f", sapply(resultados, function(r) r$stress)), ")")

paneles <- lapply(resultados, function(r) r$plot + guides(color = "none"))

## Disposicion: las cuatro primeras metricas en una cuadricula 2x2 a la
## izquierda y Jensen-Shannon a la derecha, como en la version previa.
cuadricula <- plot_grid(plotlist = paneles[1:4], labels = etiquetas[1:4],
                        label_size = 10, hjust = 0, label_x = 0.01, ncol = 2)
panel_e <- plot_grid(paneles[[5]], labels = etiquetas[5],
                     label_size = 10, hjust = 0, label_x = 0.01)
## E se deja del mismo alto que un panel de la cuadricula y centrado
## verticalmente, para no estirar su escala respecto a los demas.
col_e <- plot_grid(NULL, panel_e, NULL, ncol = 1, rel_heights = c(0.5, 1, 0.5))
grid_paneles <- plot_grid(cuadricula, col_e, ncol = 2, rel_widths = c(2, 1.05))

leyenda <- get_plot_component(resultados[[1]]$plot, "guide-box-bottom")

figura <- plot_grid(grid_paneles, leyenda, ncol = 1, rel_heights = c(1, 0.07)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("BetaDiversity.png", plot = figura, path = outpath_analisis,
       width = 30, height = 16, dpi = 300, units = "cm", bg = "white")
ggsave("BetaDiversity.png", plot = figura, path = outpath_tesis,
       width = 30, height = 16, dpi = 300, units = "cm", bg = "white")

cat("Listo 4.3\n")
