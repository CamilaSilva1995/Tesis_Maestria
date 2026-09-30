# Primeros datos de Solena: chile, maíz y tomate

Primer conjunto de datos entregado por Solena, anterior a los de fresa. Son muestras de suelo de
cultivos de chile, maíz y tomate, clasificadas con Kraken y Bracken, con las que se buscaba
detectar *Clavibacter michiganensis*. Se analizan en `Analisis_Comparativo/Clavi_Solena/`.

| Archivo o subcarpeta | Contenido |
|---|---|
| `reports/` | Reportes de Kraken por muestra. El nombre incluye el cultivo (`CHI`, `MAI`, `TOM`). |
| `kraken_braken/` | Salidas de Bracken con las abundancias reestimadas por especie. |
| `biom/` | Archivos BIOM por cultivo, construidos con los scripts `chile.sh`, `maiz.sh` y `tomate.sh`. |
| `data_reports.biom` | BIOM conjunto con todas las muestras. |
| `fastp_metadat.csv` | Estadísticas de calidad por muestra generadas con fastp. |
| `fastp_metadata_Readme.txt` | Significado de cada columna del archivo anterior (lecturas antes y después del recorte, Q30, lecturas de baja calidad, con N o demasiado cortas, duplicación). Aplica también a los `fastp_kraken_summary.csv` de `../fresa_solena/`. |
| `ALSG93AGenus.csv` | Tabla de ejemplo del libro *Statistical Analysis of Microbiome Data with R* (Xia et al.), usada para probar el análisis de potencia. |
