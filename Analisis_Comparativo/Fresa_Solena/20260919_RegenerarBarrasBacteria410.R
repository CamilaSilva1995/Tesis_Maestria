## Figura 4.10 - Abundancia relativa de bacterias por nivel taxonomico (Img/cap3/Barras_Bacteria10.png)
## Equivale a Abundance_barras(..., 10.0)[[2]] de 20230227_Funciones&Graficas.R, con las muestras
## separadas por tratamiento en paneles y la paleta comun de _paleta_tesis.R.
source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_paleta_tesis.R")
source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_barras_por_nivel.R")
datos <- cargar_fresa()
bac <- subset_taxa(datos$fil, Kingdom == "Bacteria")
cat("Bacterias:", ntaxa(bac), "taxones\n")
fig <- figura_barras_niveles(bac, umbral = 10)
guardar_figura("Barras_Bacteria10.png", fig, out_cap3, ancho = 30, alto = 20)
