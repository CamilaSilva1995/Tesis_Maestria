# Comparación Taxonomica

Queremos poder comparar a cualquier nivel taxonomico el resultado de la clasificación taxonomica de distintos clasificadores.

Queremos optimizar el ensamblaje metagenómico y las métricas de clasificación taxonómica. Por lo tanto, creamos un algoritmo de comparación entre los resultados de diferentes clasificadores taxonómicos como Kraken.

Hay varios clasificadores taxonómicos de secuencias genómicas (como Kraken) donde comparamos lo bien que clasifican a un nivel taxonómico deseado. 

Se utilizan genomas reales y datos metagenómicos simulados para evaluar la precisión de los clasificadores y ensambladores taxonómicos.


Creamos un script, con un input de datos de clasificadores metagenómicos, como por ejemplo kraken y kaiju, este genera un dataframe con la columna principal de reads y una columna con su ID taxonómico para cada uno de los clasificadores que se ingresen, con el fin de comparar el resultado entre los mismos.









