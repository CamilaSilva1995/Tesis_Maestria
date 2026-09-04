# ============================================================================
# run_all.R
# Orquesta el pipeline completo en orden. Ejecutar desde la raíz del proyecto:
#   source("Calculos_R/run_all.R")
# ============================================================================

## Resolver la raíz del proyecto de forma relativa
if (!exists("PROJECT_ROOT")) {
  wd <- getwd()
  PROJECT_ROOT <- if (basename(wd) == "Calculos_R") dirname(wd) else wd
}

scripts <- c(
  "00_setup.R",
  "01_preprocesamiento.R",
  "02_diversidad_alfa.R",
  "03_diversidad_beta.R",
  "04_exploracion_taxonomica.R",
  "05_pruebas_hipotesis.R"
)

for (s in scripts) {
  message("\n===== Ejecutando ", s, " =====")
  source(file.path(PROJECT_ROOT, "Calculos_R", s))
}

## Registrar el entorno exacto de ejecución (reproducibilidad)
writeLines(capture.output(sessionInfo()),
           file.path(PROJECT_ROOT, "Calculos_R", "Output", "sessionInfo.txt"))

message("\nPipeline completo. Figuras en Calculos_R/Results_img/, ",
        "tablas en Calculos_R/Output/.")
