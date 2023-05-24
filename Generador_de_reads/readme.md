# PyMetaSeem - Generador de reads
PyMetaSeem es un algoritmo creado desde cero que simula datos metagenómicos a partir de datos genómicos. Este simulador de datos metagenómicos fue creado a partir de reads recortadas de datos genómicos. Los genomas reales se utilizan para calificar la precisión de los clasificadores y ensambladores taxonómicos. Está pensado como un entorno conda, escalable, reproducible, fácil de usar y gratuito para el público ya que otros simuladores conocidos (como CAMISIM) tienen dificultades de instalación e implementación.




#### Crear un archivo de genomas como CAMISIM, dado el documento 'ejemplo.fasta'
1. lectura de archivo texto
 * Leer el .fasta de los contigs del genoma 
 * Ponerlo en una lista
 * Convertirlo en diccionario
 * Cortar en varios reads (por cada archivo se toman varios cortes, de contigs alazar.
 * Mezclar los cortes
 * 
2. escritura de archivos
 * Escribir en forma .fasta
 * crear el archivo .fasta con reads del metagenoma,cada uno on un indicador 
3. visualizacion grafica del resultado

queremos que el script nos arroje una lista de genomas que incluyan una lista de reads.

Esto con el fin de comprender mas a fondo el funcionamiento de CAMISIM, y tambien poder generar datos de ejemplo para probar la medida N50.
