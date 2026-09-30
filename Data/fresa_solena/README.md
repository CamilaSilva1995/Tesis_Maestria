# Datos de rizósfera de fresa (Solena)

Metagenomas *shotgun* de rizósfera de fresa (*Fragaria × ananassa*) proporcionados por la empresa
Solena Ag, con Obed Ramírez Sánchez como contacto. La clasificación taxonómica con Kraken fue
realizada por Solena; aquí se guardan los reportes por muestra, las estadísticas de calidad de
fastp y las tablas BIOM construidas con `kraken-biom`.

| Conjunto | Muestras | Categorías | Archivo BIOM |
|---|---|---|---|
| [`Data1/`](Data1/) | 58 originales, 53 tras el filtro de calidad | `healthy` (35) y `wilted` (18) | `fresa_kraken.biom`, y `fresa_genus.biom` aglomerado a género |
| [`Data2/`](Data2/) | 27 | `crop` (10), `uncultivatedLand` (10), `forestDegraded` (7) | `fresa_kraken.biom` |
| [`Data3/`](Data3/) | 8 | `crop` (4), `native` (4) | `crop_vs_native_kraken.biom` |
| [`Data_all/`](Data_all/) | Data1 + Data2 | cinco categorías: las de Data2 más `healthy` y `wilted` de Data1 | `fresa_kraken_all.biom` |

`Data1` es el conjunto que se analiza en la tesis y en el artículo. Su README documenta la descarga
desde Drive, la copia al servidor y la construcción del BIOM. Los otros tres conjuntos se
exploraron en los reportes 08, 09 y en `20230502_NuevosDatos.R` de `Analisis_Comparativo/Fresa_Solena/`.

## Filtro de calidad

Se eliminan las muestras con menos de 25 millones de lecturas después del filtrado de fastp. En
`Data1` esto excluye a MP2079, MP2080, MP2088, MP2109 y MP2137. La columna de referencia es
`Reads_A` de `fastp_kraken_summary.csv`; el significado de cada columna está en
`../solena/fastp_metadata_Readme.txt`.
