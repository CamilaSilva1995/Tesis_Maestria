# Salidas de ejemplo de Kraken y Kaiju

Fragmentos de las salidas de dos clasificadores taxonómicos sobre las mismas lecturas, usados para
desarrollar el comparador de `Comparacion_taxonomica/`.

| Archivo | Contenido |
|---|---|
| `10lineas_kraken.tsv` | Primeras diez líneas de la salida de Kraken: estado de clasificación (C/U), identificador de la lectura, taxón asignado, longitud y mapeo de k-meros. |
| `10lineas_kaiju.tsv` | Primeras diez líneas de la salida de Kaiju: estado, identificador de la lectura y taxón asignado. |
| `10lineascola_kraken.tsv`, `10lineascola_kaiju.tsv` | Últimas diez líneas de cada salida. |

Los dos programas devuelven columnas distintas; el comparador las unifica en una tabla con una fila
por lectura y una columna con el identificador taxonómico asignado por cada clasificador.
