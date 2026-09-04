# ============================================================================
# 03_diversidad_beta.R
# Diversidad beta: cinco métricas de distancia/disimilitud sobre abundancias
# relativas, visualizadas con NMDS (equivalente a la Figura 2 del artículo).
# Métricas: Bray-Curtis, Jaccard, Euclidiana, Manhattan y Jensen-Shannon (JSD).
# ============================================================================

if (!exists("percentages_fil")) source(file.path(PROJECT_ROOT, "Calculos_R", "01_preprocesamiento.R"))

## Vector de distancias a evaluar (nombre interno -> etiqueta legible)
distancias <- c(
  bray      = "Bray-Curtis",
  jaccard   = "Jaccard",
  euclidean = "Euclidiana",
  manhattan = "Manhattan",
  jsd       = "Jensen-Shannon"
)

plots <- list()
for (metodo in names(distancias)) {
  ord <- ordinate(physeq = percentages_fil, method = "NMDS", distance = metodo)
  plots[[metodo]] <- plot_ordination(
    physeq = percentages_fil, ordination = ord, color = "Treatment"
  ) +
    ggtitle(distancias[[metodo]]) +
    theme_bw()
}

## Composición de las cinco ordenaciones en una sola figura
fig_beta <- (plots[["bray"]]      | plots[["jaccard"]]) /
            (plots[["euclidean"]] | plots[["manhattan"]]) /
            (plots[["jsd"]])
fig_beta
save_fig("BetaDiversity.png", plot = fig_beta, height = 24)

message("03_diversidad_beta.R completado.")
