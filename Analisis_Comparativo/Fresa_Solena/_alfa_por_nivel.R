## Funcion comun para las Figuras 4.8 (eucariotas) y 4.11 (bacterias): diversidad alfa por nivel
## taxonomico (equivale a Alpha_diversity() de 20230227_Funciones&Graficas.R, con plot_richness),
## dibujada como caja + puntos con la paleta comun. Imprime medias por grupo y nivel.
suppressMessages(library("reshape2"))
figura_alfa_niveles <- function(ps) {
  niveles <- c(Phylum = "A) Filo", Family = "B) Familia", Genus = "C) Género", Species = "D) Especie")
  panel <- function(tax) {
    glom <- tax_glom(ps, taxrank = tax)
    alfa <- estimate_richness(glom, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
    rownames(alfa) <- sample_names(glom)          # estimate_richness altera los nombres con guion
    alfa$Tratamiento <- tratamiento_de(glom)
    larga <- melt(alfa, id.vars = "Tratamiento", measure.vars = c("Observed", "Chao1", "Shannon", "Simpson"),
                  variable.name = "Indice", value.name = "Valor")
    levels(larga$Indice) <- c("Observado", "Chao1", "Shannon", "Simpson")
    cat(niveles[tax], "- taxones:", ntaxa(glom), "\n")
    m <- aggregate(. ~ Tratamiento, alfa[, c("Observed", "Chao1", "Shannon", "Simpson", "Tratamiento")], mean)
    m[, -1] <- round(m[, -1], 3); print(m)
    ggplot(larga, aes(x = Tratamiento, y = Valor, colour = Tratamiento)) +
      geom_boxplot(outlier.shape = NA, width = 0.5) +
      geom_jitter(width = 0.12, height = 0, size = 1.3, alpha = 0.7) +
      facet_wrap(~ Indice, scales = "free_y", nrow = 1) +
      scale_colour_manual(values = colores_trat) +
      labs(x = NULL, y = NULL, colour = "Tratamiento", title = niveles[tax]) +
      tema_tesis(11) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8), axis.text.y = element_text(size = 7),
            plot.margin = margin(4, 6, 4, 6))
  }
  paneles <- lapply(names(niveles), panel)
  leyenda <- get_plot_component(paneles[[1]], "guide-box-bottom")
  grid <- plot_grid(plotlist = lapply(paneles, function(p) p + theme(legend.position = "none")), ncol = 2)
  ylab <- ggdraw() + draw_label("Medida de diversidad alfa", angle = 90, size = 12)
  plot_grid(plot_grid(ylab, grid, ncol = 2, rel_widths = c(0.03, 1)), leyenda, ncol = 1, rel_heights = c(1, 0.06))
}
