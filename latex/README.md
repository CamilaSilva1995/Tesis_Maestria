# Tesis en LaTeX

Texto de la tesis. El documento principal es `main.tex`, que incluye los capítulos de `Capitulos/`
y los apéndices de `Apendices/`. La salida es `main.pdf`.

## Estructura

| Archivo | Contenido |
|---|---|
| `main.tex` | Portada, agradecimientos, resumen, introducción, e inclusión de capítulos, apéndices y bibliografía. |
| `Capitulos/cap1_Intro.tex` | Capítulo 1, Preliminares: microbioma, metagenómica, fresa, clasificadores, formatos, índices de diversidad, pruebas de hipótesis, normalización y redes bayesianas. |
| `Capitulos/cap2.tex` | Capítulo 2, Objetivos. |
| `Capitulos/cap3.tex` | Capítulo 3, Metodología del análisis exploratorio. |
| `Capitulos/cap4.tex` | Capítulo 4, Resultados. |
| `Capitulos/cap5.tex` | Capítulo 5, Discusión y conclusiones. |
| `Capitulos/cap6.tex` | Capítulo 6, Perspectivas y trabajo a futuro: redes bayesianas y simulador. |
| `Apendices/Ap.tex`, `Bp.tex` | Apéndices de códigos e imágenes. Por completar. |
| `ref.bib` | Bibliografía en BibTeX, estilo `apalike`. |
| `Img/` | Figuras por capítulo (`cap1/` a `cap4/`) y logotipos de la portada. |

Los archivos `.aux`, `.log`, `.toc`, `.out`, `.bbl` y `.synctex.gz` son productos de la
compilación. Los `.bak` son copias de versiones anteriores.

## Compilar

```bash
cd latex
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

Se necesitan los paquetes `babel` (español), `natbib`, `hyperref`, `graphicx`, `listings`,
`svg` y `amsmath`. El `.gitignore` de esta carpeta excluye la carpeta `.oberon/`.

## Figuras

Las figuras del capítulo 4 están en `Img/cap3/` y se regeneran a 300 dpi con los scripts
`20260919_Regenerar*.R` y `20260924_Regenerar*.R` de `Analisis_Comparativo/Fresa_Solena/`. La
tabla de correspondencia entre script y figura está en el README de esa carpeta.

## Ecuaciones

Todas las ecuaciones en display están numeradas por capítulo y tienen etiqueta para citarlas con
`\ref`. Las principales son `eq:chao1`, `eq:shannon`, `eq:simpson`, `eq:bray-curtis`,
`eq:jaccard-cuant`, `eq:jensen-shannon`, `eq:t-statistic`, `eq:welch-statistic`,
`eq:estadistico-f`, `eq:mann-whitney`, `eq:matriz-conteos` y `eq:simplex`.
