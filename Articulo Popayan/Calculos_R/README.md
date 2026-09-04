# Cálculos en R — Análisis estadístico de la diversidad microbiana

Este directorio contiene el **flujo computacional reproducible** que soporta el artículo
*"Análisis estadístico de la diversidad microbiana a distintos niveles taxonómicos en
datos metagenómicos de alta dimensión"*.

El objetivo es que cualquier persona pueda **reproducir los cálculos e imágenes** del
artículo a partir de la matriz de conteos taxonómicos (formato BIOM) generada con Kraken
y los metadatos de las muestras, siguiendo un pipeline modular y documentado.

## Diseño del pipeline

El análisis se estructura como una secuencia de scripts numerados. Cada script es un
**módulo independiente** que carga su entrada, ejecuta una etapa del análisis y persiste
sus salidas (figuras en `Results_img/`, tablas intermedias en `Output/`). Esta separación
por etapas facilita la depuración, la trazabilidad y la reejecución selectiva.

```
Data/
  fresa_kraken.biom      # matriz de conteos taxón x muestra (salida de kraken-biom)
  metadata.csv           # metadatos de las muestras (incluye la variable Treatment)

Calculos_R/
  00_setup.R                 # paquetes, rutas y utilidades comunes
  01_preprocesamiento.R      # carga BIOM -> objeto phyloseq, limpieza y control de calidad
  02_diversidad_alfa.R       # índices Chao1, Shannon y Simpson
  03_diversidad_beta.R       # distancias (Bray-Curtis, Jaccard, Euclidiana, Manhattan, JSD) + NMDS
  04_exploracion_taxonomica.R# aglomeración por nivel taxonómico y barras de abundancia
  05_pruebas_hipotesis.R     # t de Student (var. iguales/distintas) y Mann-Whitney/Wilcoxon
  run_all.R                  # ejecuta todo el pipeline en orden
```

## Requisitos

- **R** >= 4.2
- Paquetes de Bioconductor: `phyloseq`
- Paquetes de CRAN: `ggplot2`, `vegan`, `dplyr`, `patchwork`, `RColorBrewer`, `stringi`

Instalación (una sola vez):

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("phyloseq")
install.packages(c("ggplot2", "vegan", "dplyr", "patchwork", "RColorBrewer", "stringi"))
```

## Cómo generar la matriz BIOM (etapa previa, en Bash)

Los datos crudos fueron clasificados taxonómicamente con **Kraken** y reestimados con
**Bracken**. La matriz de conteos se consolida en formato BIOM con `kraken-biom`:

```bash
conda activate metagenomics
kraken-biom kraken_results/* --fmt json -o Data/fresa_kraken.biom
```

## Cómo ejecutar

Desde la raíz del proyecto (o desde este directorio), en R:

```r
source("Calculos_R/run_all.R")
```

O paso a paso, en orden numérico (`00_setup.R` primero). Todas las rutas se resuelven de
forma relativa a través de `00_setup.R`, por lo que no es necesario editar rutas absolutas.

## Notas de reproducibilidad

- La semilla aleatoria se fija en `00_setup.R` (`set.seed(1995)`) para que las ordenaciones
  NMDS y cualquier submuestreo sean reproducibles.
- El filtro de control de calidad elimina 5 muestras con < 25 millones de lecturas tras el
  filtrado (`MP2079`, `MP2080`, `MP2088`, `MP2109`, `MP2137`), dejando **53 muestras**
  (35 saludables y 18 no saludables).
- `sessionInfo()` se guarda en `Output/sessionInfo.txt` al final de `run_all.R` para dejar
  registro exacto de las versiones utilizadas.

> Este flujo es una versión curada y depurada de los scripts de investigación originales
> ubicados en `../Fresa_Solena/`, reorganizados para publicación y reproducibilidad.
