# Generalizando el N50: métricas para ensamblajes de metagenomas

El N50 es la medida más usada para evaluar la contigüidad de un ensamblaje genómico. Los contigs se
ordenan de mayor a menor y se suman sus longitudes; el N50 es la longitud del contig con el que la
suma acumulada alcanza el 50 % del tamaño total del ensamblaje. El N90 se define igual con el 90 %.

En metagenómica un ensamblaje mezcla contigs de muchos genomas de tamaños muy distintos, de modo
que un solo N50 puede quedar dominado por los organismos más abundantes o de genoma más grande.
Esta carpeta evalúa si el N50 sigue siendo pertinente para metagenomas y explora cómo generalizarlo.

| Archivo | Contenido |
|---|---|
| `N50.sh` | Script de Bash que calcula estadísticas básicas de un ensamblaje en FASTA: número de contigs, tamaño total, contig más largo, N50 y N90. Uso: `bash N50.sh ensamblaje.fasta`. |
| `N50.ipynb` | Implementación en Python del N50 genómico y exploración de su generalización a metagenomas. |

## Exploración de la generalización

El cuaderno parte de un archivo de ejemplo con siete lecturas, lo separa en dos genomas con
lecturas tomadas al azar y compara el N50 del metagenoma completo con el N50 de cada genoma. El
resultado muestra el problema: el N50 del conjunto coincide con el de uno de los genomas y no
refleja al otro.

Las definiciones alternativas que se plantean son:

1. Promediar los N50 de cada genoma presente.
2. Calcular la lista de N50 por genoma y tomar el N50 de esa lista.

Ambas requieren saber a qué genoma pertenece cada contig, es decir, una clasificación taxonómica
previa, y estimar el tamaño del genoma de cada taxón. Para probar las definiciones con una
composición conocida se usan los metagenomas simulados con PyMetaSeem (`Generador_de_reads/`).

Esta línea quedó en fase exploratoria y no forma parte de los resultados de la tesis.
