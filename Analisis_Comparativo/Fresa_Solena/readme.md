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
| `20230314_Redes.R` | Redes simples con phyloseq y redes de coocurrencia con *Fusarium* como género de interés. Los umbrales 0.8 y 0.5 que usa producen el grafo completo; ver `20261002_RegenerarRedesMuestras32.R`. Los tres scripts `20261002_RegenerarComplementarios3*.R` usan el cargador común `_cargar_complementarios.R`, que une los conjuntos Data_all y Data3 con sus metadatos por identificador; reemplazan a `20230320`, `20230321`, `20230327` y `20230502`. |
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

Cada script reproduce una figura a 300 dpi con rótulos en español y la guarda en `latex/Img/` (carpeta `cap1`, `cap2` o `cap3` según el capítulo)
y en `Results_img/`. Todas usan la misma paleta de `_paleta_tesis.R`. El número al final del nombre de cada script es el de la figura en la tesis. Todos los scripts con
componente aleatoria usan `set.seed(2026)`, la semilla única de la tesis.

| Script | Figura | Archivo |
|---|---|---|
| `20260919_RegenerarPruebatExplicacion13.R` | 1.3 | `pruebat_explicacion.png` (en `latex/Img/cap1/`) |
| `20261002_RegenerarRedesMuestras32.R` | 3.2 | `Redes_muestras.png` (en `latex/Img/cap2/`) |
| `20261002_RegenerarRedesCoocurrencia33.R` | 3.3 y Tabla 3.1 | `Redes_coocurrencia.png` (en `latex/Img/cap2/`), `Redes_metricas.csv` |
| `20261002_RegenerarComplementarios34.R` | 3.4 | `Complementarios_alfa.png` (en `latex/Img/cap2/`), `Complementarios_kruskal.csv` |
| `20261002_RegenerarComplementarios35.R` | 3.5 y Tabla 3.2 | `Complementarios_beta.png` (en `latex/Img/cap2/`), `Complementarios_permanova*.csv` |
| `20261002_RegenerarComplementarios36.R` | 3.6 | `Complementarios_nativo.png` (en `latex/Img/cap2/`) |
| `20260924_RegenerarAlphaCrudosFiltrados41.R` | 4.1 y Tabla 4.1 | `AlphaDiversity_CrudosFiltrados.png` |
| `20260924_RegenerarBarrasLecturas42.R` | 4.2 | `Barras.png` |
| `20260924_RegenerarBetaDiversity43.R` | 4.3 | `BetaDiversity.png` |
| `20260930_RegenerarPermanovaPermdisp44.R` | 4.4 y Tabla 4.2 | `PERMANOVA_PERMDISP.png`, `PERMANOVA_PERMDISP_resultados.csv` |
| `20260919_RegenerarBarrasFilo45.R` | 4.5 | `BarrasFilo.png` |
| `20260906_RegenerarBarrasFiloBN46.R` | 4.6 | `BarrasFilo_blackandwhite.png` |
| `20260924_RegenerarBarrasEukarya47.R` | 4.7 | `Barras_Eukarya10.png` |
| `20260919_RegenerarAlphaEukarya48.R` | 4.8 | `Alpha_Eukarya.png` |
| `20260919_RegenerarBetaEukarya49.R` | 4.9 | `Beta_Eukarya.png`, `Beta_Eukarya_ordinaciones.rds` |
| `20260919_RegenerarBarrasBacteria410.R` | 4.10 | `Barras_Bacteria10.png` |
| `20260919_RegenerarAlphaBacteria411.R` | 4.11 | `Alpha_Bacteria.png` |
| `20260919_RegenerarBetaBacteria412.R` | 4.12 | `Beta_Bacteria.png`, `Beta_Bacteria_ordinaciones.rds` |
| `20260919_RegenerarTTestShannon413.R` | 4.13 y Tabla 4.2 | `tTest_Shannon.png` |
| `20260919_RegenerarTTestFusarium414.R` | 4.14 y Tabla 4.2 | `tTest_Shannon_Fusarium.png` |
| `20260919_RegenerarWilcoxonShannon415.R` | 4.15 y Tabla 4.2 | `Wilcoxon_Shannon.png` |
| `20261001_RegenerarRarefaccionNormalizacion416.R` | 4.16 y Tablas 4.3, 4.4 | `Rarefaccion_Normalizacion.png`, `Rarefaccion_pruebas.csv`, `Normalizacion_PERMANOVA.csv` |
| `20261001_RegenerarProporcionesFusarium417.R` | 4.17 y Tabla 4.5 | `Proporciones_Fusarium.png`, `Proporciones_generos_candidatos.csv` |

Archivos auxiliares que los scripts cargan con `source()`:

| Archivo | Función |
|---|---|
| `_paleta_tesis.R` | Paleta única de la tesis (saludable `#F8766D`, no saludable `#00BFC4`, un color fijo por filo, paleta cualitativa para familias, géneros y especies), tema de ggplot2, cargador del conjunto Data1 con el filtro de calidad y función para guardar a 300 dpi. |
| `_cargar_complementarios.R` | Cargador de los conjuntos Data_all y Data3 con sus metadatos y colores de las cinco categorías. |
| `_barras_por_nivel.R`, `_alfa_por_nivel.R`, `_beta_por_nivel.R` | Funciones comunes de las figuras por nivel taxonómico (4.7 a 4.12). |
| `_prueba_t_figura.R` | Figura común de las pruebas t (4.13 y 4.14): histogramas por grupo y distribuciones t con la región de rechazo. |

El Apéndice A de la tesis (`latex/Apendices/Ap.tex`) reproduce esta correspondencia figura → script → datos.

`20261001_RegenerarRarefaccionNormalizacion416.R` calcula las curvas de rarefacción, rarifica a la profundidad mínima con semilla 2026, repite las pruebas sobre los índices rarificados y compara el PERMANOVA con abundancias relativas, TMM (edgeR) y datos rarificados; guarda `Results_img/Rarefaccion_pruebas.csv` y `Normalizacion_PERMANOVA.csv`. `20260930_RegenerarPermanovaPermdisp44.R` aplica PERMANOVA (`adonis2`, Bray-Curtis y Jaccard) y PERMDISP (`betadisper`) sobre la composición completa, guarda la tabla `Results_img/PERMANOVA_PERMDISP_resultados.csv` y la figura con la ordenación PCoA y las distancias al centroide.

## Otros archivos

- `Results_img/`: todas las figuras generadas, incluidas las versiones antiguas y las de los datos adicionales.
- `percentages_df.csv`: tabla de abundancias relativas en formato largo, exportada desde phyloseq.
- `Presentacion.qmd` y `Presentacion1.qmd`: presentación en Quarto del análisis. La segunda se presentó en el Congreso Colombiano de Matemáticas 2025.

## Resultados principales

- El filtro de calidad es indispensable: sin él, las muestras con muy pocas lecturas ocultan cualquier diferencia.
- Las plantas saludables muestran mayor diversidad alfa en todos los índices y niveles, pero la diferencia es pequeña y no es significativa en la comunidad completa (Welch p = 0.1501, Mann-Whitney p = 0.2586).
- La diversidad beta no separa los grupos con ninguna medida.
- La diversidad de Shannon de los géneros de eucariotas, grupo que incluye a *Fusarium*, tampoco difiere entre grupos (medias 3.077 y 3.035; t con varianzas iguales p = 0.260, Welch p = 0.378, Mann-Whitney p = 0.963). El valor p = 0.0017 que daba el script de 2023 era un artefacto: unía el índice con los metadatos por posición después de `psmelt`/`reshape`, que reordenan las muestras. Script 414.
- Las plantas no saludables son más heterogéneas entre sí: la prueba F rechaza la igualdad de varianzas de Shannon (p < 0.001).
- Conjuntos complementarios: el tipo de sitio explica el 33 % de la variación en composición (PERMANOVA p = 0.001); cultivo saludable y marchito no difieren entre sí pero ambos difieren de los suelos cercanos, no cultivados y de bosque degradado (R² de 0.22 a 0.33). Scripts 34 a 36.
- Las redes de muestras (Jaccard) no separan los grupos (fracción de aristas dentro del grupo 52-56 %, esperada 54 %, p > 0.3). En la red SparCC entre géneros, *Fusarium* no tiene vecinos con |r| ≥ 0.7 y sus vecinos más fuertes (Alnitak) son hongos de su propio linaje.
- La proporción de *Fusarium* (0.14 % de las lecturas en ambos grupos, p = 0.98) y la de los géneros candidatos no difieren entre grupos; la única diferencia nominal, *Phytophthora* (p = 0.013, mayor en saludables), no sobrevive a la corrección de Benjamini-Hochberg. Script 417, tabla `Results_img/Proporciones_generos_candidatos.csv`.
- Tras rarificar a profundidad común, Shannon no cambia (p = 0.150) y la diferencia de riqueza entre grupos se reduce a menos de la mitad: parte de la mayor riqueza de las saludables era efecto de la profundidad.
- El PERMANOVA no detecta diferencia de composición promedio (R² = 0.024, p = 0.17), pero PERMDISP confirma que las plantas no saludables son más dispersas en composición (p = 0.023).
