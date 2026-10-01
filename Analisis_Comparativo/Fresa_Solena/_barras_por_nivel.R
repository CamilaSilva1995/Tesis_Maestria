## Funcion comun para las Figuras 4.7 (eucariotas) y 4.10 (bacterias): abundancia relativa por
## nivel taxonomico con aglomeracion al 10 % (los taxones que no superan ese umbral en ninguna
## muestra se agrupan en "Otros"), muestras separadas por tratamiento y paleta comun.
figura_barras_niveles <- function(ps, umbral = 10) {
  panel <- function(tax, titulo) {
    df <- psmelt(transform_sample_counts(tax_glom(ps, taxrank = tax), function(x) 100 * x / sum(x)))
    df$Taxon <- df[[tax]]
    visibles <- unique(df$Taxon[df$Abundance >= umbral])
    df$Taxon <- ifelse(df$Taxon %in% visibles, df$Taxon, "Otros")
    df <- aggregate(Abundance ~ Sample + Tratamiento + Taxon, df, sum)
    orden <- names(sort(tapply(df$Abundance, df$Taxon, mean), decreasing = TRUE))
    orden <- c(setdiff(orden, "Otros"), "Otros")
    df$Taxon <- factor(df$Taxon, levels = orden)
    df$Sample <- factor(df$Sample, levels = sort(unique(df$Sample)))
    ggplot(df, aes(x = Sample, y = Abundance, fill = Taxon)) +
      geom_col(width = 0.9) +
      facet_grid(. ~ Tratamiento, scales = "free_x", space = "free_x") +
      scale_fill_manual(values = paleta_taxones(orden)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
      labs(x = NULL, y = NULL, fill = NULL, title = titulo) +
      tema_tesis(11) +
      theme(axis.text.x = element_text(angle = 90, size = 4.5, hjust = 1, vjust = 0.5),
            legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 8),
            plot.margin = margin(4, 6, 4, 6)) +
      guides(fill = guide_legend(nrow = ceiling(length(orden) / 4)))
  }
  paneles <- list(panel("Phylum", "A) Filo"), panel("Family", "B) Familia"),
                  panel("Genus", "C) Género"), panel("Species", "D) Especie"))
  grid <- plot_grid(plotlist = paneles, ncol = 2)
  ylab <- ggdraw() + draw_label("Abundancia relativa (%)", angle = 90, size = 12)
  xlab <- ggdraw() + draw_label("Muestra", size = 12)
  plot_grid(plot_grid(ylab, grid, ncol = 2, rel_widths = c(0.03, 1)), xlab, ncol = 1, rel_heights = c(1, 0.04))
}
