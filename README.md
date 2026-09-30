# Tesis de Maestría en Ciencias Matemáticas

Repositorio de la tesis *Análisis estadístico de la diversidad microbiana a distintos niveles taxonómicos en el microbioma rizosférico de plantas de fresa saludables y no saludables*, del Posgrado Conjunto en Ciencias Matemáticas UMSNH-UNAM (Morelia, México).
Autora: Paula Camila Silva Gómez. Asesora: Dra. Nelly Sélem Mojica.

El trabajo compara el microbioma de la rizósfera de plantas de fresa saludables y no saludables a
partir de metagenomas *shotgun* proporcionados por la empresa Solena Ag, usando índices de
diversidad alfa y beta, exploración taxonómica por niveles y pruebas de hipótesis. Como líneas
complementarias incluye herramientas para comparar clasificadores taxonómicos, una generalización
de la métrica N50 para metagenomas y un simulador de lecturas metagenómicas (PyMetaSeem).

## Mapa del repositorio

| Carpeta | Contenido |
|---|---|
| [`latex/`](latex/) | Texto de la tesis en LaTeX, figuras y bibliografía. Se compila a `latex/main.pdf`. |
| [`Analisis_Comparativo/`](Analisis_Comparativo/) | Análisis exploratorio y estadístico en R. La subcarpeta `Fresa_Solena/` contiene el análisis principal de la tesis. |
| [`Data/`](Data/) | Todos los datos: tablas de abundancia (BIOM), metadatos, reportes de Kraken, genomas y lecturas. |
| [`Articulo_Popayan/`](Articulo_Popayan/) | Artículo derivado de la tesis para la Revista Científica (Universidad Distrital), con figuras, cálculos y documentos de envío. |
| [`Comparacion_taxonomica/`](Comparacion_taxonomica/) | Scripts en Python para comparar las salidas de distintos clasificadores taxonómicos (Kraken, Kaiju). |
| [`N50/`](N50/) | Cálculo de la métrica N50 para ensamblajes y exploración de su generalización a metagenomas. |
| [`Generador_de_reads/`](Generador_de_reads/) | PyMetaSeem, simulador de lecturas metagenómicas a partir de genomas de referencia. |
| [`Redes/`](Redes/) | Espacio reservado para el modelado con redes bayesianas, trabajo en curso. |

## Flujo principal de la tesis

1. **Datos.** Los reportes de clasificación con Kraken se convierten a BIOM con `kraken-biom`
   (ver [`Data/fresa_solena/Data1/Readme.md`](Data/fresa_solena/Data1/Readme.md)).
2. **Análisis.** Los cuadernos R Markdown numerados en
   [`Analisis_Comparativo/Fresa_Solena/`](Analisis_Comparativo/Fresa_Solena/) cargan el BIOM como
   objeto phyloseq, filtran por calidad, calculan diversidades y aplican las pruebas de hipótesis.
3. **Figuras.** Los scripts `20260919_Regenerar*.R` y `20260924_Regenerar*.R` de esa misma carpeta
   producen las figuras finales a 300 dpi y las guardan en `latex/Img/cap3/`.
4. **Tesis.** Se compila con `pdflatex` y `bibtex` desde `latex/` (ver [`latex/README.md`](latex/README.md)).

## Requisitos

- R (4.x) con los paquetes `phyloseq`, `vegan`, `ggplot2`, `dplyr`, `plyr`, `patchwork` y `edgeR`.
- Python 3 con `pandas`, `numpy` y `ete3` para las herramientas de comparación taxonómica y el simulador.
- Un ambiente conda con `kraken2` y `kraken-biom` para regenerar las tablas BIOM.
- LaTeX (`pdflatex`, `bibtex`) para la tesis y `pandoc` con LibreOffice para el artículo.

Los datos ocupan cerca de 1.3 GB y están versionados en el repositorio.
