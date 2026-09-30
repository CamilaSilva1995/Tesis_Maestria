# Data2: muestras por tipo de sitio

Segundo conjunto de datos entregado por Solena. Son 27 muestras de suelo clasificadas según el
sitio de muestreo: `crop` (10, cultivo de fresa), `uncultivatedLand` (10, tierra sin cultivar) y
`forestDegraded` (7, bosque degradado).

| Archivo | Contenido |
|---|---|
| `fresa_kraken.biom` | Tabla de abundancias construida con `kraken-biom` a partir de los reportes de Kraken. |
| `metadata.csv` | Identificador de muestra y categoría. |
| `fastp_kraken_summary.csv` | Estadísticas de calidad por muestra (lecturas antes y después del filtrado, Q30, duplicación, porcentaje clasificado). |
| `taxonomic_profiling/` | Reportes de Kraken por muestra. |

Se exploró en `Analisis_Comparativo/Fresa_Solena/08_ExploracionDatosNuevos.Rmd` y, unido a
`Data1`, en `09_ExploracionDatosTotal.Rmd`. Ver `../Data_all/`.
