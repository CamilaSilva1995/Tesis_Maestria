## Funcion comun para las Figuras 4.13 y 4.14: histograma del indice de Shannon por grupo y
## distribuciones t con la region de rechazo (alfa = 0.05, bilateral) y el estadistico observado,
## para la prueba con varianzas iguales y para la de Welch. Mismo estilo que la Figura 1.3.
figura_prueba_t <- function(shannon, tratamiento, titulo_x = "Índice de diversidad de Shannon", bins = 25) {
  df <- data.frame(Valor = shannon, Tratamiento = tratamiento)
  medias <- aggregate(Valor ~ Tratamiento, df, mean)
  cat("Medias por grupo:\n"); print(medias)
  cat("Desviaciones estandar:\n"); print(aggregate(Valor ~ Tratamiento, df, sd))
  pA <- ggplot(df, aes(x = Valor, fill = Tratamiento)) +
    geom_histogram(bins = bins, colour = "white", linewidth = 0.2) +
    geom_vline(data = medias, aes(xintercept = Valor), linetype = "dashed", linewidth = 0.7, colour = "grey20") +
    facet_grid(Tratamiento ~ .) +
    scale_fill_manual(values = colores_trat, guide = "none") +
    labs(x = titulo_x, y = "Número de muestras", title = "A) Distribución por grupo (línea: media)") +
    tema_tesis(11)
  x <- split(df$Valor, df$Tratamiento)
  t1 <- t.test(x[[1]], x[[2]], var.equal = TRUE); t2 <- t.test(x[[1]], x[[2]], var.equal = FALSE)
  print(t1); print(t2)
  panel_t <- function(tt, titulo) {
    gl <- as.numeric(tt$parameter); est <- as.numeric(tt$statistic); crit <- qt(0.975, gl)
    cur <- data.frame(x = seq(-4.5, 4.5, length = 600)); cur$y <- dt(cur$x, gl)
    ggplot(cur, aes(x, y)) +
      geom_area(data = subset(cur, x >= -crit & x <= crit), fill = "#C6E2FF", alpha = 0.9) +
      geom_area(data = subset(cur, x < -crit), fill = "grey70") +
      geom_area(data = subset(cur, x > crit), fill = "grey70") +
      geom_line(colour = "grey40", linewidth = 0.4) +
      geom_vline(xintercept = est, colour = "royalblue4", linewidth = 1) +
      annotate("text", x = est, y = max(cur$y) * 0.55, label = sprintf("t = %.3f", est), colour = "royalblue4",
               angle = 90, vjust = -0.5, size = 3.6) +
      annotate("text", x = c(-crit, crit), y = -0.02, label = sprintf("%.3f", c(-crit, crit)), size = 3) +
      labs(title = titulo, x = sprintf("Distribución t con %.1f grados de libertad; p = %.3f", gl, tt$p.value),
           y = "Densidad") +
      tema_tesis(11) + theme(panel.grid = element_blank())
  }
  pB <- panel_t(t1, "B) Varianzas iguales; región de rechazo en gris (α = 0.05)")
  pC <- panel_t(t2, "C) Varianzas distintas (Welch)")
  list(figura = plot_grid(pA, plot_grid(pB, pC, ncol = 1), ncol = 2, rel_widths = c(1, 1.05)), t_igual = t1, t_welch = t2)
}
