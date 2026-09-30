# Data

Todos los datos usados por las herramientas y análisis del repositorio. Cada subcarpeta tiene su
propio README con el origen de los datos y los comandos con que se procesaron. En total ocupan
cerca de 1.3 GB y están versionados en git.

| Subcarpeta | Tamaño | Contenido | Usado en |
|---|---|---|---|
| [`fresa_solena/`](fresa_solena/) | 349 MB | Datos principales de la tesis: metagenomas *shotgun* de rizósfera de fresa proporcionados por Solena. Cuatro conjuntos: `Data1` (saludable / no saludable), `Data2` (tres tipos de sitio), `Data3` (cultivo / nativo) y `Data_all` (unión con cinco categorías). | `Analisis_Comparativo/Fresa_Solena/`, `Articulo_Popayan/` |
| [`fresa_paper/`](fresa_paper/) | 440 MB | Lecturas públicas del NCBI del estudio de Yang et al. (2020) sobre fresa con y sin oídio, con sus reportes de Kraken y el BIOM. | `Analisis_Comparativo/Fresa_paper/` |
| [`solena/`](solena/) | 60 MB | Primeros datos de Solena: reportes de Kraken y Bracken de cultivos de chile, maíz y tomate, y sus BIOM. | `Analisis_Comparativo/Clavi_Solena/` |
| [`Clavibacter/`](Clavibacter/) | 234 MB | Genomas de *Clavibacter* completos y recortados, para probar el simulador de lecturas y la métrica N50. | `Generador_de_reads/`, `N50/` |
| [`ReadsFonty/`](ReadsFonty/) | 185 MB | Diez genomas de referencia de una comunidad simulada (bacterias y levaduras), con sus versiones recortadas y perfiles de abundancia. | `Generador_de_reads/` |
| [`kraken_kaiju/`](kraken_kaiju/) | 56 KB | Fragmentos de las salidas de Kraken y Kaiju para desarrollar el comparador de clasificadores. | `Comparacion_taxonomica/` |

## Convenciones

- Las tablas de abundancia se guardan en formato BIOM (JSON) generado con `kraken-biom` y se leen
  en R con `phyloseq::import_biom`.
- Cada conjunto de fresa tiene un `metadata.csv` con el identificador de muestra y su categoría, y
  un `fastp_kraken_summary.csv` con las estadísticas de calidad por muestra usadas para el filtro.
- Los identificadores de muestra de Solena tienen la forma `MP2079`, `MD2055`, etc.
