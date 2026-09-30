# Clavibacter en cultivos de chile, maíz y tomate (Solena)

Primer ejercicio de análisis con datos de Solena, anterior al trabajo con fresa. Se buscaba
detectar la presencia de la bacteria fitopatógena *Clavibacter michiganensis* en muestras de
suelo de tres cultivos y comparar las comunidades microbianas entre ellos.

Los datos de entrada son los reportes de Kraken de `Data/solena/reports/` y el archivo
`Data/solena/data_reports.biom`.

| Archivo | Contenido |
|---|---|
| `chile.sh`, `maiz.sh`, `tomate.sh` | Construyen un archivo BIOM por cultivo con `kraken-biom` a partir de los reportes de cada muestra. |
| `trim-clavi.sh` | Extrae de un reporte de Kraken las líneas correspondientes a *Clavibacter michiganensis* y las guarda en `trim-reports/`. Uso: `bash trim-clavi.sh muestra.report`. |
| `data_reports.R` | Carga el BIOM conjunto como objeto phyloseq, siguiendo el tutorial de metagenómica de The Carpentries, y explora las diversidades. |
| `alpha_beta_diversity(phyloseq).R` | Diversidad alfa y beta por cultivo calculadas con phyloseq. |
| `alpha_beta_diversity(vegan).R` | Las mismas diversidades calculadas con vegan, para comparar ambas implementaciones. |
| `alphadiversity_*_plot.pdf`, `betadiversity_*_plot.pdf` | Gráficas resultantes para chile, maíz y tomate. |
| `paneles.png` | Panel resumen con las gráficas de los tres cultivos. |

Este análisis sirvió para aprender el flujo Kraken → BIOM → phyloseq que después se aplicó a los
datos de fresa. No forma parte de los resultados de la tesis.
