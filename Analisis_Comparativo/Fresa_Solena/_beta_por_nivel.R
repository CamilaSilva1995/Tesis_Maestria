## Funcion comun para las Figuras 4.9 (eucariotas) y 4.12 (bacterias): NMDS con Bray-Curtis por
## nivel taxonomico. Equivale a Beta_diversity(..., 'bray') de 20230227_Funciones&Graficas.R, sin
## las etiquetas de muestra (ilegibles con 53 puntos) y con el estres de cada ordenacion en el
## titulo. Semilla 2026 antes de cada NMDS; las ordenaciones se guardan en un .rds para poder
## reajustar la figura sin recalcular (borrar el archivo fuerza el recalculo).
suppressMessages(library("vegan"))
figura_beta_niveles <- function(ps, cache) {
  niveles <- c(Phylum = "A) Filo", Family = "B) Familia", Genus = "C) Género", Species = "D) Especie")
  if (file.exists(cache)) { ords <- readRDS(cache); cat("Ordenaciones leidas de", cache, "\n") } else {
    ords <- lapply(names(niveles), function(tax) {
      rel <- transform_sample_counts(tax_glom(ps, taxrank = tax), function(x) 100 * x / sum(x))
      set.seed(2026); o <- ordinate(rel, method = "NMDS", distance = "bray", trymax = 50, trace = 0)
      list(scores = data.frame(scores(o, display = "sites")), stress = o$stress, muestras = sample_names(rel))
    }); names(ords) <- names(niveles); saveRDS(ords, cache)
  }
  trat <- tratamiento_de(ps)
  panel <- function(tax) {
    o <- ords[[tax]]; sc <- o$scores; sc$Tratamiento <- trat[o$muestras]
    cat(niveles[tax], "- estrés:", round(o$stress, 3), if (o$stress < 0.01) "(ordenación degenerada: muy pocos taxones)" else "", "\n")
    ggplot(sc, aes(NMDS1, NMDS2, colour = Tratamiento)) +
      geom_point(size = 2.2, alpha = 0.85) +
      { if (o$stress > 0.01) stat_ellipse(level = 0.68, linewidth = 0.5, show.legend = FALSE) } +
      scale_colour_manual(values = colores_trat) +
      labs(title = sprintf("%s (estrés = %.3f)", niveles[tax], o$stress), colour = "Tratamiento") +
      tema_tesis(11) + theme(axis.text = element_text(size = 8))
  }
  paneles <- lapply(names(niveles), panel)
  leyenda <- get_plot_component(paneles[[1]], "guide-box-bottom")
  grid <- plot_grid(plotlist = lapply(paneles, function(p) p + theme(legend.position = "none")), ncol = 2)
  plot_grid(grid, leyenda, ncol = 1, rel_heights = c(1, 0.06))
}
