## Figura 4.12 - Diversidad beta (NMDS, Bray-Curtis) de bacterias por nivel taxonomico
##                (Img/cap3/Beta_Bacteria.png). Hasta 2026 esta figura provenia de
##                20230227_Funciones&Graficas.R a 960x540 px.
source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_paleta_tesis.R")
source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_beta_por_nivel.R")
datos <- cargar_fresa()
bac <- subset_taxa(datos$fil, Kingdom == "Bacteria")
fig <- figura_beta_niveles(bac, file.path(out_res, "Beta_Bacteria_ordinaciones.rds"))
guardar_figura("Beta_Bacteria.png", fig, out_cap3, ancho = 30, alto = 19)
