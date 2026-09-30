# Comparación de clasificadores taxonómicos

Herramientas en Python para comparar, a cualquier nivel taxonómico, las asignaciones que distintos
clasificadores de secuencias (Kraken, Kaiju) hacen sobre las mismas lecturas. La idea es evaluar
qué tan bien clasifica cada programa usando genomas reales y metagenomas simulados con PyMetaSeem,
en los que la respuesta correcta se conoce de antemano.

Cada clasificador devuelve una tabla con columnas distintas. El comparador las unifica en una sola
tabla con una fila por lectura y una columna con el identificador taxonómico asignado por cada
programa, para poder contar coincidencias y discrepancias por nivel.

| Archivo | Contenido |
|---|---|
| `Proyecto1.ipynb` | Desarrollo paso a paso del comparador con pandas: lectura de las salidas de Kraken y Kaiju, eliminación de las columnas que no son comunes y unión por identificador de lectura. |
| `Proyecto_metanomica.ipynb`, `Proyecto_metanomica.py` | Versión ordenada del comparador, con detección automática del tipo de salida (Kraken o Kaiju). |
| `get_full_lineages.ipynb`, `get_full_lineages.py` | Dada una lista de identificadores taxonómicos, obtiene el linaje completo (reino a especie) consultando la taxonomía del NCBI con `ete3`, y lo guarda en `full_lineages.tsv`. Necesario para comparar a un nivel taxonómico distinto del asignado. |
| `full_lineages.tsv` | Linajes obtenidos para los taxones de las salidas de ejemplo. |
| `cutout.ipynb`, `cutout.py` | Función `cutout(read, i, n)` que recorta una lectura desde la posición `i` con longitud `n`. Se usa para simplificar encabezados FASTA y para generar lecturas cortas de prueba. |

Los datos de ejemplo están en `Data/kraken_kaiju/`. Para `get_full_lineages.py` se requiere
`ete3` con la base de datos de taxonomía del NCBI descargada (`NCBITaxa().update_taxonomy_database()`).

Esta línea de trabajo se planteó como parte de la optimización de métricas de clasificación
taxonómica y quedó en fase exploratoria; la tesis se centró en el análisis de la fresa.
