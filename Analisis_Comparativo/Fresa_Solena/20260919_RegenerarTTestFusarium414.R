## Figura 4.14 - Prueba t de Student sobre el indice de Shannon de los generos de eucariotas
##                (Img/cap3/tTest_Shannon_Fusarium.png). Reproduce 20230410_PruebasdeHipotesisMedias.R
##                (lineas 282-390): subconjunto Eukaryota aglomerado a genero (58 generos, entre
##                ellos Fusarium), Shannon con vegan::diversity, prueba t y prueba de Welch.
## CORRECCION (2026-09-29): el script de 2023 unia el indice con los metadatos por posicion despues
## de psmelt/reshape, que reordenan las muestras; cada valor quedaba asignado al tratamiento de otra
## muestra y daba p = 0.0017. Aqui el indice se calcula directamente sobre la tabla de conteos del
## objeto aglomerado, cuyas columnas conservan el orden de sam_data, y se verifica con stopifnot.
source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_paleta_tesis.R")
source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_prueba_t_figura.R")
suppressMessages(library("vegan"))
datos <- cargar_fresa()
glom <- tax_glom(subset_taxa(datos$fil, Kingdom == "Eukaryota"), taxrank = "Genus")
cat("Generos de eucariotas:", ntaxa(glom), "\n")
otu <- t(glom@otu_table@.Data)
shannon <- diversity(otu, "shannon")
trat <- tratamiento_de(glom)
stopifnot(identical(names(shannon), names(trat)))
res <- figura_prueba_t(shannon, trat, titulo_x = "Índice de Shannon (géneros de eucariotas)")
cat("Mann-Whitney:\n"); print(wilcox.test(shannon[trat == "Saludable"], shannon[trat == "No saludable"]))
guardar_figura("tTest_Shannon_Fusarium.png", res$figura, out_cap3, ancho = 30, alto = 16)
