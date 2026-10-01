## Figura 4.7 - Abundancia relativa de eucariotas por nivel taxonomico (Img/cap3/Barras_Eukarya10.png)
## Equivale a Abundance_barras(..., 10.0)[[2]] de 20230227_Funciones&Graficas.R, con las muestras
## separadas por tratamiento en paneles y la paleta comun de _paleta_tesis.R.
source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_paleta_tesis.R")
source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_barras_por_nivel.R")
datos <- cargar_fresa()
euk <- subset_taxa(datos$fil, Kingdom == "Eukaryota")
cat("Eucariotas:", ntaxa(euk), "taxones\n")
fig <- figura_barras_niveles(euk, umbral = 10)
guardar_figura("Barras_Eukarya10.png", fig, out_cap3, ancho = 30, alto = 20)
