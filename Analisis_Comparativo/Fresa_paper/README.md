# Réplica con los datos públicos de Yang et al. (2020)

Análisis de los datos públicos del artículo
[Comparison of the rhizosphere soil microbial community structure and diversity between powdery mildew-infected and noninfected strawberry plants in a greenhouse by high-throughput sequencing technology](https://link.springer.com/article/10.1007/s00284-020-01948-x)
(*Current Microbiology*, 2020), que compara la rizósfera de fresa con y sin oídio. Se usó como
referencia metodológica antes de analizar los datos de Solena.

Los datos de entrada están en `Data/fresa_paper/`: cuatro muestras descargadas del NCBI con SRA
Toolkit, clasificadas con Kraken 2 y convertidas a `fresa_paper.biom`. Son datos de amplicón 16S,
por lo que la clasificación con Kraken es solo aproximada; el flujo recomendado para estos datos es
QIIME 2.

| Reporte | Contenido |
|---|---|
| `230227_Reporte1Exploracion.Rmd` | Preprocesamiento, filtro de calidad y diversidades alfa y beta sobre el conjunto completo y el filtrado. |
| `230306_Reporte2Funciones.Rmd` | Subconjuntos de Eukaryota y Bacteria por nivel taxonómico, con las funciones automatizadas. |
| `230313_Reporte3PruebasdeHipotesis.Rmd` | Pruebas de hipótesis sobre los índices de diversidad entre plantas infectadas y no infectadas. |

Con solo cuatro muestras los resultados no son concluyentes. El valor de este ejercicio fue
establecer el flujo de trabajo y las funciones que después se reutilizaron en `Fresa_Solena/`.
