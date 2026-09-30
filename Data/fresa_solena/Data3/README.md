# Data3: cultivo frente a suelo nativo

Tercer conjunto de datos entregado por Solena, con ocho muestras: cuatro de cultivo de fresa
(`crop`) y cuatro de suelo nativo (`native`).

| Archivo | Contenido |
|---|---|
| `crop_vs_native_kraken.biom` | Tabla de abundancias construida con `kraken-biom`. |
| `crop_vs_native_rawCounts_S.csv`, `.tsv`, `.xlsx` | La misma tabla de conteos crudos a nivel de especie, en formato plano, tal como la entregó Solena. |
| `crop_vs_native_rawCounts_S_metadata.csv`, `.xlsx` | Metadatos de las muestras entregados por Solena. |
| `metadata.csv` | Identificador de muestra y categoría, en el formato que usan los scripts de R. |
| `fastp_kraken_summary.csv` | Estadísticas de calidad por muestra. |
| `kraken_results/` | Reportes de Kraken por muestra. |
| `Data3-20230503T004053Z-001.zip` | Paquete original descargado de Drive. |

Se exploró en `Analisis_Comparativo/Fresa_Solena/20230502_NuevosDatos.R`. La tabla `.csv` es la
que carga el cuaderno `20230428_ML.ipynb`.
