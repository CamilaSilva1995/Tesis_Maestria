## Figura 4.13 - Prueba t de Student sobre el indice de Shannon de la comunidad completa
##                (Img/cap3/tTest_Shannon.png). Reproduce 20230410_PruebasdeHipotesisMedias.R
##                (lineas 30-116): Shannon con vegan::diversity sobre la tabla de conteos de las
##                53 muestras filtradas, prueba t con varianzas iguales y prueba de Welch.
##                Hasta 2026 la figura era un montaje manual de tres imagenes a 960x540 px.
source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_paleta_tesis.R")
source("/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/_prueba_t_figura.R")
suppressMessages(library("vegan"))
datos <- cargar_fresa()
otu <- t(datos$fil@otu_table@.Data)
shannon <- diversity(otu, "shannon")
trat <- tratamiento_de(datos$fil)
stopifnot(identical(names(shannon), names(trat)))
res <- figura_prueba_t(shannon, trat)
cat("Shapiro-Wilk por grupo:\n"); print(tapply(shannon, trat, function(v) shapiro.test(v)$p.value))
cat("Prueba F de varianzas:\n"); print(var.test(shannon[trat == "Saludable"], shannon[trat == "No saludable"]))
cat("Mann-Whitney:\n"); print(wilcox.test(shannon[trat == "Saludable"], shannon[trat == "No saludable"]))
guardar_figura("tTest_Shannon.png", res$figura, out_cap3, ancho = 30, alto = 16)
