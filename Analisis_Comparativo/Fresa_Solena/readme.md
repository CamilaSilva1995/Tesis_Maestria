# Análisis exploratorio de datos metagenómicos de fresa (Solena)

Objetivo: encontrar características diferenciadoras entre los microbiomas rizosféricos de plantas
de fresa saludables y no saludables, a partir de metagenomas *shotgun* proporcionados por
[Solena](https://solena.ag). Tras el filtro de calidad quedan 53 muestras: 35 saludables
(`healthy`) y 18 no saludables (`wilted`).

El análisis se divide en tres partes: exploración con diversidades alfa y beta y composición por
nivel taxonómico; validación con pruebas de hipótesis; y una exploración preliminar de redes de
coocurrencia y de clasificación con aprendizaje automático.

**Datos de entrada:** `Data/fresa_solena/Data1/fresa_kraken.biom` y `metadata.csv`. Los datos
adicionales por tipo de cultivo están en `Data2/`, `Data3/` y `Data_all/`.

## Reportes (R Markdown)

Los cuadernos numerados son la versión ordenada y comentada del análisis. Varios tienen su PDF compilado junto al `.Rmd`.

| Reporte | Contenido |
|---|---|
| `01_Exploracion.Rmd` | Preprocesamiento, descripción del formato BIOM y del objeto phyloseq, filtro de calidad (25 millones de lecturas), diversidad alfa y beta con todas las medidas sobre el conjunto completo, y una primera visualización con redes simples. No hay separación clara entre grupos. |
| `02_ExploracionSubconjuntos.Rmd` | Subconjuntos por reino (Bacteria y Eukaryota) y por nivel taxonómico, con sus diversidades alfa y beta. |
| `03_FuncionesAutomatizacion.Rmd` | Funciones que automatizan la aglomeración por nivel taxonómico y las gráficas de barras, alfa y beta. |
| `04_PruebasdeHipotesis.Rmd` | Pruebas t de Student (varianzas iguales y de Welch), F sobre varianzas y Wilcoxon-Mann-Whitney sobre Shannon y Chao1, en el conjunto completo y en subconjuntos. |
| `05_FusariumActinobacteria.Rmd` | Barras de abundancia y diversidad beta para *Fusarium* (género) y Actinobacteria (filo) como taxones de interés. |
| `06_DiversidadesAlfa10%.Rmd` | Diversidad alfa con aglomeración al 10 % para mejorar la visualización. |
| `07_RedesCoocurrencia.Rmd` | Redes de coocurrencia con MicNet (UMAP, SparCC) y Alnitak, tomando *Fusarium* como taxón principal. Las salidas están en `Data/fresa_solena/Data1/Redes/`. |
| `08_ExploracionDatosNuevos.Rmd` | Exploración del segundo conjunto de datos, con tres categorías por tipo de cultivo. |
| `09_ExploracionDatosTotal.Rmd` | Exploración de los dos conjuntos unidos, con cinco categorías. |
| `10_Normalización.Rmd` | Normalización de la tabla de conteos con edgeR. |
| `11_ExploraciondDatosNormalizados.Rmd` | Repetición del análisis principal sobre los datos normalizados. |
| `12_RarefacciónSaturacionDeMuestra.Rmd` | Curvas de rarefacción: las muestras están saturadas, hay lecturas suficientes. |
| `13_FuncionesAutomatizacionNormalizados.Rmd` | Funciones automatizadas aplicadas a los datos normalizados. |

## Scripts de desarrollo (2023)

Scripts de R con fecha, en el orden en que se hizo el análisis. Los reportes anteriores los recopilan.

| Script | Contenido |
|---|---|
| `20230130_PreprocesamientoDatos.R` | Descarga de los datos, creación del BIOM y carga como objeto phyloseq. Primera vista de las diversidades. |
| `20230213_DiversidadesAlfa&Beta.R` | Barras de abundancia, porcentajes y diversidades por nivel taxonómico. |
| `20230220_Funciones.R` | Tres funciones para aglomerar por nivel y graficar barras, alfa y beta. |
| `20230227_Funciones&Graficas.R` | Versión mejorada de las funciones; genera todas las gráficas por filo, familia, género y especie, separando Bacteria y Eukaryota. |
| `20230306_Potencia&PruebaHipotesis.R` | Análisis de potencia siguiendo a Xia et al., *Statistical Analysis of Microbiome Data with R*. |
| `20230314_Redes.R` | Redes simples con phyloseq y redes de coocurrencia con *Fusarium* como género de interés. |
| `20230320_NuevosDatos.R` | Segundo conjunto de datos, tres categorías por sitio de muestreo. |
| `20230321_NuevosDatosAll.R` | Unión de los dos conjuntos: cinco categorías por tipo de cultivo y estado. |
| `20230327_PruebasOrden.R` | Reordenamiento de las cinco categorías para las gráficas de barras. |
| `20230403_Actinobacteria.R` | Barras y diversidad beta para Actinobacteria. |
| `20230403_Oomycota&Fusarium.R` | Barras y diversidad beta para Oomycota y *Fusarium*. |
| `20230410_PruebasdeHipotesisMedias.R` | Prueba t sobre las medias del índice de Shannon y prueba de normalidad de Shapiro-Wilk. |
| `20230411_Preprocesamiento&Normalización.R` | Normalización con edgeR. |
| `20230419_Rarefaccion.R` | Curvas de rarefacción como medida de saturación de las muestras. |
| `20230425_Rarefaccion(Normalizacion).R` | Rarefacción como método de normalización. |
| `20230427_PruebasdeHipotesisVarianzas.R` | Prueba F sobre las varianzas de Chao1. |
| `20230428_ML.ipynb` | Primer intento de clasificación saludable frente a no saludable con scikit-learn. Exploratorio. |
| `20230502_NuevosDatos.R` | Tercer conjunto de datos: cultivo frente a nativo. |
| `20230510_PruebaWilcoxon.R` | Prueba de Wilcoxon-Mann-Whitney sobre Shannon (W = 376, p = 0.2586). |
| `Mann-Whitney.R` | Versión resumida de la prueba anterior. |

## Scripts que regeneran las figuras de la tesis (2026)

Cada script reproduce una figura a 300 dpi con rótulos en español y la guarda en `latex/Img/cap3/`
y en `Results_img/`. El número al final del nombre es la figura de la tesis.

| Script | Figura | Archivo |
|---|---|---|
| `20260919_RegenerarPruebatExplicacion13.R` | 1.3 | `pruebat_explicacion.png` |
| `20260924_RegenerarAlphaCrudosFiltrados41.R` | 4.1 | `AlphaDiversity_CrudosFiltrados.png` |
| `20260924_RegenerarBarrasLecturas42.R` | 4.2 | `Barras.png` |
| `20260924_RegenerarBetaDiversity43.R` | 4.3 | `BetaDiversity.png` |
| `20260919_RegenerarBarrasFilo44.R` | 4.4 | `BarrasFilo.png` |
| `20260924_RegenerarBarrasEukarya46.R` | 4.6 | `Barras_Eukarya10.png` |
| `20260919_RegenerarBarrasBacteria47.R` | 4.7 | `Barras_Bacteria10.png` |
| `20260919_RegenerarAlphaEukarya48.R` | 4.8 | `Alpha_Eukarya.png` |
| `20260919_RegenerarAlphaBacteria49.R` | 4.9 | `Alpha_Bacteria.png` |
| `20260919_RegenerarTTestFusarium413.R` | 4.13 | `tTest_Shannon_Fusarium.png` |
| `20260919_RegenerarWilcoxonShannon414.R` | 4.14 | `Wilcoxon_Shannon.png` |

`20260906_RegenerarBarrasFilo.R` es una versión anterior del script de la figura 4.4.

## Otros archivos

- `Results_img/`: todas las figuras generadas, incluidas las versiones antiguas y las de los datos adicionales.
- `percentages_df.csv`: tabla de abundancias relativas en formato largo, exportada desde phyloseq.
- `Presentacion.qmd` y `Presentacion1.qmd`: presentación en Quarto del análisis. La segunda se presentó en el Congreso Colombiano de Matemáticas 2025.

## Resultados principales

- El filtro de calidad es indispensable: sin él, las muestras con muy pocas lecturas ocultan cualquier diferencia.
- Las plantas saludables muestran mayor diversidad alfa en todos los índices y niveles, pero la diferencia es pequeña y no es significativa en la comunidad completa (Welch p = 0.1501, Mann-Whitney p = 0.2586).
- La diversidad beta no separa los grupos con ninguna medida.
- La diversidad de Shannon de los géneros de eucariotas, grupo que incluye a *Fusarium*, sí difiere significativamente (p = 0.0017 con varianzas iguales, p = 0.0130 con Welch).
- Las plantas no saludables son más heterogéneas entre sí: la prueba F rechaza la igualdad de varianzas de Shannon (p < 0.001).
