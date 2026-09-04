# ============================================================================
# 01_preprocesamiento.R
# Carga la matriz BIOM como objeto phyloseq, limpia la taxonomía, incorpora
# metadatos y aplica el control de calidad. Produce dos objetos:
#   - fresa_kraken      : datos crudos (58 muestras)
#   - fresa_kraken_fil  : datos filtrados por calidad (53 muestras)
# ============================================================================

if (!exists("PROJECT_ROOT")) source(file.path(getwd(), "Calculos_R", "00_setup.R"))

## --- Carga de la matriz de conteos (BIOM) ---------------------------------
# El BIOM es una matriz dispersa (taxón x muestra) con metadatos asociados.
fresa_kraken <- import_biom(BIOM_FILE)

## --- Limpieza de la tabla de taxonomía ------------------------------------
# 1) nombrar las columnas con los rangos taxonómicos
colnames(fresa_kraken@tax_table@.Data) <- TAX_RANKS
# 2) quitar el prefijo de rango (p. ej. "p__Bacteria" -> "Bacteria")
fresa_kraken@tax_table@.Data <- substr(fresa_kraken@tax_table@.Data, 4, 100)

## --- Homologar identificadores de muestra ---------------------------------
# recortar los nombres de columna de la tabla de OTU para que coincidan con
# los identificadores de los metadatos (6 caracteres, p. ej. "MD2055")
colnames(fresa_kraken@otu_table@.Data) <- substr(colnames(fresa_kraken@otu_table@.Data), 1, 6)

## --- Metadatos ------------------------------------------------------------
metadata_fresa <- read.csv2(METADATA_FILE, header = FALSE, row.names = 1, sep = ",")
fresa_kraken@sam_data <- sample_data(metadata_fresa)
fresa_kraken@sam_data$Sample <- row.names(fresa_kraken@sam_data)
colnames(fresa_kraken@sam_data) <- c("Treatment", "Samples")

## --- Control de calidad ---------------------------------------------------
# Se descartan las muestras con < 25 M lecturas tras el filtrado de calidad.
fresa_kraken_fil <- prune_samples(
  !(sample_names(fresa_kraken) %in% SAMPLES_TO_REMOVE),
  fresa_kraken
)

## --- Verificación ---------------------------------------------------------
message("Muestras antes del filtro:  ", nsamples(fresa_kraken))       # 58
message("Muestras después del filtro: ", nsamples(fresa_kraken_fil))  # 53
message("Distribución por tratamiento (filtradas):")
print(table(as.data.frame(fresa_kraken_fil@sam_data)$Treatment))       # 35 healthy / 18 wilted

## --- Composición por reino (tamaño del espacio taxonómico) ----------------
n_bacteria <- sum(fresa_kraken_fil@tax_table@.Data[, "Kingdom"] == "Bacteria")
n_eukarya  <- sum(fresa_kraken_fil@tax_table@.Data[, "Kingdom"] == "Eukaryota")
message("OTU Bacteria: ", n_bacteria, " | OTU Eukaryota: ", n_eukarya)

## --- Abundancias relativas (normalización composicional) ------------------
# Cada muestra se lleva al símplex (sus proporciones suman 100).
percentages     <- transform_sample_counts(fresa_kraken,     function(x) x * 100 / sum(x))
percentages_fil <- transform_sample_counts(fresa_kraken_fil, function(x) x * 100 / sum(x))

## Persistir tabla larga de porcentajes para etapas posteriores
percentages_df <- psmelt(percentages_fil)
write.csv(percentages_df, file.path(OUTPUT_DIR, "percentages_df.csv"), row.names = FALSE)

message("01_preprocesamiento.R completado.")
