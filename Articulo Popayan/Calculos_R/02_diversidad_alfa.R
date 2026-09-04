# ============================================================================
# 02_diversidad_alfa.R
# Índices de diversidad alfa (Observado, Chao1, Shannon, Simpson) sobre los
# datos crudos y filtrados, comparando muestras saludables y no saludables.
# Requiere que 01_preprocesamiento.R se haya ejecutado antes.
# ============================================================================

if (!exists("fresa_kraken_fil")) source(file.path(PROJECT_ROOT, "Calculos_R", "01_preprocesamiento.R"))

## --- Tabla de índices de diversidad alfa ----------------------------------
index <- estimate_richness(
  fresa_kraken_fil,
  measures = c("Observed", "Chao1", "Shannon", "Simpson")
)
write.csv(index, file.path(OUTPUT_DIR, "alpha_diversity_index.csv"))

## --- Visualización: datos crudos vs. filtrados ----------------------------
p_crudos <- plot_richness(
  physeq = fresa_kraken,
  measures = c("Observed", "Chao1", "Shannon", "Simpson"),
  x = "Treatment", color = "Treatment"
) + ggtitle("A) Datos crudos")

p_filtrados <- plot_richness(
  physeq = fresa_kraken_fil,
  measures = c("Observed", "Chao1", "Shannon", "Simpson"),
  x = "Treatment", color = "Treatment"
) + ggtitle("B) Datos filtrados por calidad")

## Figura combinada (equivalente a la Figura 1 del artículo)
(p_crudos / p_filtrados)
save_fig("AlphaDiversity_CrudosFiltrados.png")

message("02_diversidad_alfa.R completado.")
