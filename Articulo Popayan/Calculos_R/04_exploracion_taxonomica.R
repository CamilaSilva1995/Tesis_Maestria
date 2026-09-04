# ============================================================================
# 04_exploracion_taxonomica.R
# Exploración multinivel: aglomeración por nivel taxonómico y barras de
# abundancia (absoluta y relativa), separando por reino cuando corresponde.
# Genera insumos equivalentes a las Figuras 3 y 4 del artículo.
# ============================================================================

if (!exists("fresa_kraken_fil")) source(file.path(PROJECT_ROOT, "Calculos_R", "01_preprocesamiento.R"))

## --- Función de aglomeración + barras a un nivel taxonómico dado ----------
# rank      : uno de TAX_RANKS (p. ej. "Phylum", "Genus")
# top_frac  : conserva los taxones cuya abundancia relativa supera este umbral
#             (agrupando el resto como "Others") para mejorar la legibilidad.
barras_por_nivel <- function(physeq, rank, top_frac = 0.10, etiqueta = rank) {
  glom <- tax_glom(physeq, taxrank = rank, NArm = FALSE)
  rel  <- transform_sample_counts(glom, function(x) x * 100 / sum(x))

  # colapsar taxones poco abundantes en "Others"
  df <- psmelt(rel)
  df[[rank]] <- as.character(df[[rank]])
  df[[rank]][df$Abundance < top_frac * 100] <- "Others"

  ggplot(df, aes(x = Sample, y = Abundance, fill = .data[[rank]])) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Treatment, scales = "free_x", space = "free_x") +
    ylab("Abundancia relativa (%)") + xlab("Muestra") +
    ggtitle(etiqueta) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "bottom")
}

## --- Barras a nivel de filo (todas las muestras) --------------------------
fig_filo <- barras_por_nivel(fresa_kraken_fil, "Phylum", top_frac = 0.0, etiqueta = "Filo")
save_fig("BarrasFilo.png", plot = fig_filo)

## --- Subconjuntos por reino ----------------------------------------------
bacteria <- subset_taxa(fresa_kraken_fil, Kingdom == "Bacteria")
eukarya  <- subset_taxa(fresa_kraken_fil, Kingdom == "Eukaryota")

## --- Composición bacteriana multinivel (Figura 4) -------------------------
fig_b_phylum  <- barras_por_nivel(bacteria, "Phylum",  etiqueta = "A) Filo")
fig_b_family  <- barras_por_nivel(bacteria, "Family",  etiqueta = "B) Familia")
fig_b_genus   <- barras_por_nivel(bacteria, "Genus",   etiqueta = "C) Género")
fig_b_species <- barras_por_nivel(bacteria, "Species", etiqueta = "D) Especie")

fig_bacteria <- (fig_b_phylum | fig_b_family) / (fig_b_genus | fig_b_species)
save_fig("Barras_Bacteria10.png", plot = fig_bacteria, height = 24)

## --- Composición eucariota multinivel (opcional) --------------------------
fig_e_phylum  <- barras_por_nivel(eukarya, "Phylum",  etiqueta = "A) Filo")
fig_e_family  <- barras_por_nivel(eukarya, "Family",  etiqueta = "B) Familia")
fig_e_genus   <- barras_por_nivel(eukarya, "Genus",   etiqueta = "C) Género")
fig_e_species <- barras_por_nivel(eukarya, "Species", etiqueta = "D) Especie")

fig_eukarya <- (fig_e_phylum | fig_e_family) / (fig_e_genus | fig_e_species)
save_fig("Barras_Eukarya10.png", plot = fig_eukarya, height = 24)

message("04_exploracion_taxonomica.R completado.")
