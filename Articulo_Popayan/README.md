# Artículo para la Revista Científica (Universidad Distrital)

Versión en español del artículo "Análisis estadístico de la diversidad microbiana a distintos
niveles taxonómicos en datos metagenómicos de alta dimensión", ajustada a la plantilla y a las
normas de la revista (ISSN 0124-2253, E-ISSN 2344-8350) para la primera ronda de revisión.

## Archivos

| Archivo | Contenido |
|---|---|
| `Articulo_RevistaCientifica_ES.docx` | Manuscrito anonimizado listo para subir a la plataforma OJS. |
| `Articulo_RevistaCientifica_ES_vista_previa.pdf` | Vista previa generada con LibreOffice (sustituye Times New Roman por Liberation Serif). |
| `articulo_ES.md` | Fuente en Markdown del manuscrito. Editar este archivo y volver a ejecutar `build_docx.py`. |
| `build_docx.py` | Genera el .docx con pandoc y aplica el formato de la revista (Times New Roman 12, interlineado 1.5, A4, márgenes 2.5/3 cm). |
| `figuras/Figura_1.jpg` a `Figura_6.jpg` | Figuras a 300 dpi para cargar por separado en la plataforma, como exige la revista. |
| `original/Articulo_Diversidad_Microbiana_Fresa_EN_v2.docx` | Versión original en inglés recibida por el editor. |
| `original/Plantilla_RevistaCientifica.docx` | Plantilla oficial descargada de la página de la revista. |
| `original/articulo_EN.md`, `original/plantilla.md` | Texto extraído de los dos documentos anteriores. |
| `img/media/` | Imágenes extraídas de la versión original (baja resolución, solo de referencia). |

## Regenerar el .docx

```bash
python3 build_docx.py
```

## Pendientes antes de la versión final (después de la revisión)

- Restituir nombres, filiaciones, correos institucionales y ORCID de los tres autores.
- Indicar la dirección del repositorio público con el código.
- Indicar agradecimientos institucionales y fuentes de financiación.
- Completar la lista de autores de Yang et al. (2020) en formato APA 7.
