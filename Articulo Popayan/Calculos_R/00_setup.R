# ============================================================================
# 00_setup.R
# Configuración común del pipeline: paquetes, rutas relativas y utilidades.
# Se debe ejecutar (o cargar con source) antes que cualquier otro script.
# ============================================================================

## --- Paquetes -------------------------------------------------------------
suppressPackageStartupMessages({
  library("phyloseq")      # manejo de datos metagenómicos (objeto phyloseq)
  library("ggplot2")       # visualización
  library("vegan")         # índices de diversidad y distancias ecológicas
  library("dplyr")         # manipulación de data frames
  library("patchwork")     # composición de figuras
  library("RColorBrewer")  # paletas de color
  library("stringi")       # utilidades de cadenas
})

## --- Reproducibilidad -----------------------------------------------------
set.seed(1995)  # fija la semilla para NMDS y cualquier submuestreo aleatorio

## --- Rutas relativas ------------------------------------------------------
# Se resuelve la raíz del proyecto de forma robusta: este script asume que
# vive en <raiz>/Calculos_R/. Ajuste PROJECT_ROOT si cambia la estructura.
if (!exists("PROJECT_ROOT")) {
  # Cuando se ejecuta con source(), normalizePath del wd suele bastar.
  candidate <- getwd()
  if (basename(candidate) == "Calculos_R") {
    PROJECT_ROOT <- dirname(candidate)
  } else {
    PROJECT_ROOT <- candidate
  }
}

DATA_DIR   <- file.path(PROJECT_ROOT, "Data")
OUTPUT_DIR <- file.path(PROJECT_ROOT, "Calculos_R", "Output")
IMG_DIR    <- file.path(PROJECT_ROOT, "Calculos_R", "Results_img")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(IMG_DIR,    showWarnings = FALSE, recursive = TRUE)

BIOM_FILE     <- file.path(DATA_DIR, "fresa_kraken.biom")
METADATA_FILE <- file.path(DATA_DIR, "metadata.csv")

## Muestras a excluir por control de calidad (< 25 M lecturas tras filtrado)
SAMPLES_TO_REMOVE <- c("MP2079", "MP2080", "MP2088", "MP2109", "MP2137")

## Niveles taxonómicos de trabajo
TAX_RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

## --- Utilidad para guardar figuras de forma homogénea ---------------------
save_fig <- function(filename, plot = last_plot(),
                     width = 30, height = 15, dpi = 300, units = "cm") {
  ggsave(filename = filename, plot = plot, path = IMG_DIR,
         width = width, height = height, dpi = dpi, units = units)
  invisible(file.path(IMG_DIR, filename))
}

message("00_setup.R cargado. PROJECT_ROOT = ", PROJECT_ROOT)
