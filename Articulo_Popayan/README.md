# Artículo para la Revista Científica (Universidad Distrital)

Envío 25531. Versión en español del artículo "Análisis estadístico de la diversidad microbiana
a distintos niveles taxonómicos en datos metagenómicos de alta dimensión", ajustada a la
plantilla y a las normas de la revista (ISSN 0124-2253, E-ISSN 2344-8350) para la primera ronda
de revisión, según lo solicitado por el equipo editorial: Resumen, Abstract y Resumo;
Introducción; Metodología; Resultados; Conclusiones; sin nombres de autores.

## Qué subir a la plataforma OJS

1. `Articulo_RevistaCientifica_ES.docx`: manuscrito anonimizado (texto del artículo).
2. `figuras/Figura1_DiversidadAlfa.jpg` a `Figura6_Wilcoxon_Shannon.jpg`: las seis figuras a
   300 dpi, cargadas por separado como exige la revista.
3. `documentos_envio/`: carta de originalidad y formato de cesión de derechos, ya diligenciados
   y enviados con la primera versión, y la declaración de conflicto de intereses y financiación,
   que la revista exige y no se envió en la primera ronda. Completar la financiación y firmarla.

## Contenido de la carpeta

| Ruta | Contenido |
|---|---|
| `Articulo_RevistaCientifica_ES.docx` | Manuscrito listo para subir. Times New Roman 12, interlineado 1.5, A4, márgenes de la plantilla, ecuaciones numeradas, referencias APA 7 con DOI. |
| `Articulo_RevistaCientifica_ES_CON_AUTORES_no_subir.docx` | Copia de respaldo con nombres, filiaciones, correos y ORCID en notas al pie, contribuciones con nombres y repositorio. No subir en la ronda de revisión; sirve de base para la versión final. Se genera desde `articulo_ES_con_autores.md`. |
| `Articulo_RevistaCientifica_ES_vista_previa.pdf` | Vista previa generada con LibreOffice (sustituye Times New Roman por Liberation Serif). |
| `articulo_ES.md` | Fuente del manuscrito en Markdown. Editar este archivo y volver a ejecutar `build_docx.py`. |
| `build_docx.py` | Genera el .docx con pandoc y aplica el formato de la revista. |
| `figuras/` | Seis figuras en JPEG a 300 dpi, con rótulos en español. |
| `calculos/Figura1_*.R` a `Figura6_*.R` | Scripts de R que generan cada figura a partir de los datos de Solena (`Data/fresa_solena/Data1`). |
| `calculos/Verificacion_formulas.R` | Comprueba que las ecuaciones 1 a 3, las medidas de diversidad beta y los estadísticos de las pruebas coinciden con lo que calculan phyloseq y vegan. |
| `calculos/Potencia_estadistica.R` | Tamaño del efecto, potencia y PERMANOVA que respaldan la discusión sobre el tamaño de muestra. |
| `documentos_envio/` | Carta de originalidad y formato de cesión de derechos diligenciados, más la declaración de conflicto de intereses y financiación que la revista también exige, redactada y pendiente de completar y firmar. |
| `original/Articulo_Diversidad_Microbiana_Fresa_EN_v2.docx` | Versión original en inglés enviada en la primera ronda. |
| `original/envio_original_25531_submission-files.zip` | Paquete completo descargado de la plataforma con lo que se subió en la primera ronda. |
| `original/Plantilla_RevistaCientifica.docx` | Plantilla oficial descargada de la página de la revista. |
| `original/articulo_EN.md`, `original/plantilla.md`, `original/plantilla_media/` | Texto e imágenes extraídos de los dos documentos anteriores. El banner de la revista sale de aquí. |

## Regenerar el manuscrito

```bash
python3 build_docx.py                                   # version anonima
python3 build_docx.py articulo_ES_con_autores.md Articulo_RevistaCientifica_ES_CON_AUTORES_no_subir.docx
```

## Cambios respecto a la versión en inglés

- Traducción completa al español, con Resumen, Abstract y Resumo de menos de 200 palabras.
- Sin nombres, filiaciones, correos ni ORCID de los autores, por la revisión doble ciego.
- Fórmula de Simpson corregida a la forma de Gini-Simpson, que es la que calcula phyloseq.
- Chao1 descrito como estimador que usa los singletons y doubletons, no que los elimina.
- Prueba de Wilcoxon descrita como suma de rangos para muestras independientes.
- La prueba sobre *Fusarium* descrita como prueba sobre el índice de Shannon de los géneros de eucariotas. Su resultado cambió: la versión en inglés reportaba p = 0.0017 por un error de emparejamiento de muestras en el script; con las muestras bien emparejadas no hay diferencia significativa (p = 0.26). Ver `calculos/Figura5_tStudent_Fusarium.R`.
- Se agregaron los resultados de las pruebas de Shapiro-Wilk y F, que muestran que el grupo no saludable es significativamente más heterogéneo.
- Bracken retirado de la metodología, porque la matriz BIOM se construyó solo con las salidas de Kraken.
- Figuras reemplazadas por las versiones a 300 dpi con rótulos en español.
- Ecuaciones como objetos del editor de ecuaciones de Word, centradas y numeradas a la derecha.
- Texto ajustado para no superar las 20 páginas incluyendo bibliografía (límite de la revista).
- Agradecimiento a Solena Ag y a Obed Ramírez Sánchez.

## Pendientes antes de la versión final (después de la revisión)

- Restituir nombres, filiaciones, correos institucionales y ORCID de los tres autores.
- Indicar la dirección del repositorio público con el código.
- Indicar agradecimientos institucionales y fuentes de financiación.
- Completar la lista de autores de Yang et al. (2020) en formato APA 7.
