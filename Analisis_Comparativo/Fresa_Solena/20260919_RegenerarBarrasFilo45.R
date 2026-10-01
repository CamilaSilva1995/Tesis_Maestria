## Figura 4.5 - Abundancia absoluta y relativa por filo, muestras separadas por grupo
##               (Img/cap3/BarrasFilo.png)
## Reproduce el analisis de 02_ExploracionSubconjuntos.Rmd (chunks absolute_plot y relative_plot)
## con la paleta fija por filo de _paleta_tesis.R. Los 44 filos de la tabla no se pueden distinguir
## por color, asi que se muestran los 10 de mayor abundancia media y el resto se agrupa en "Otros".
##   A) lecturas asignadas a cada filo (millones)   B) porcentaje por muestra
## Las muestras aparecen agrupadas por tratamiento (paneles), ordenadas por identificador.

source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_paleta_tesis.R")
datos <- cargar_fresa()

n_top <- 10
glom_abs <- tax_glom(datos$fil, taxrank = "Phylum")
glom_rel <- tax_glom(datos$rel, taxrank = "Phylum")
df_abs <- psmelt(glom_abs); df_rel <- psmelt(glom_rel)
top <- names(sort(tapply(df_rel$Abundance, df_rel$Phylum, mean), decreasing = TRUE))[1:n_top]
cat("Filos principales (abundancia relativa media, %):\n")
print(round(sort(tapply(df_rel$Abundance, df_rel$Phylum, mean), decreasing = TRUE)[1:n_top], 2))

preparar <- function(df) {
  df$Filo <- factor(ifelse(df$Phylum %in% top, df$Phylum, "Otros"), levels = c(top, "Otros"))
  df$Sample <- factor(df$Sample, levels = sort(unique(df$Sample)))
  df
}
df_abs <- preparar(df_abs); df_rel <- preparar(df_rel)
pal <- paleta_para_filos(levels(df_abs$Filo))

barras <- function(df, y, etiqueta_y, titulo) {
  ggplot(df, aes(x = Sample, y = .data[[y]], fill = Filo)) +
    geom_col(width = 0.9) +
    facet_grid(. ~ Tratamiento, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Muestra", y = etiqueta_y, fill = "Filo", title = titulo) +
    tema_tesis() +
    theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 1, vjust = 0.5),
          legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 9))
}
df_abs$Millones <- df_abs$Abundance / 1e6
pA <- barras(df_abs, "Millones",  "Lecturas clasificadas (millones)", "A) Abundancia absoluta por filo")
pB <- barras(df_rel, "Abundance", "Abundancia relativa (%)",         "B) Abundancia relativa por filo")

leyenda <- get_plot_component(pA + guides(fill = guide_legend(nrow = 2)), "guide-box-bottom")
fig <- plot_grid(plot_grid(pA + theme(legend.position = "none"), pB + theme(legend.position = "none"), ncol = 1),
                 leyenda, ncol = 1, rel_heights = c(1, 0.08))
guardar_figura("BarrasFilo.png", fig, out_cap3, ancho = 30, alto = 24)
