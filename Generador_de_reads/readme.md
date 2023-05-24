# PyMetaSeem - Generador de reads
PyMetaSeem es un algoritmo creado desde cero que simula datos metagenómicos a partir de datos genómicos. Este simulador de datos metagenómicos fue creado a partir de reads recortadas de datos genómicos. Los genomas reales se utilizan para calificar la precisión de los clasificadores y ensambladores taxonómicos. Está pensado como un entorno conda, escalable, reproducible, fácil de usar y gratuito para el público ya que otros simuladores conocidos (como CAMISIM) tienen dificultades de instalación e implementación.


##IMAGEN WORKFLOW##

Este algoritmo está compuesto de varias funciones:

PyMetaSeem recibe un conjunto de archivos con diferentes genomas en formato FASTA, para esto se tienen las opciones de inicializar el algoritmo con numero de reads y longitud de reads ya sea manual o en un archivo de texto (txt) que contenga elnombre de los archivos y la cantidad de reads por genoma (perfil taxonomico).

Se tiene una función que se encarga de la lectura de los fasta, convirtiendolos en un diccionario que contiene cada secuencia de nucleótidos como valor y su nombre o identificador como clave, lo que facilita la gestión de los datos dentro de python. 

Teniendo los genomas en diccionarios, se crean varias funciones para tomar las longitudes de los reads a recortar y las proporciones,  tomando en cuenta el tamaño original de los contigs. Luego de obtener las longitudes y proporciones, se cortan los reads aleatorias (m) de una longitud determinada (n) a partir de la secuencia de nucleótidos.

A los reads cortados se les calcula el reverso complementario. 

Luego se tiene una funcion que llama a todas las funciones anteriores y crea un diccionario de reads y una ultima función convierte estos diccionarios en archivos FASTQ comprimidos que muestra los datos metagenómicos con su calidad respectiva.



Esto con el fin de comprender mas a fondo el funcionamiento de CAMISIM, y tambien poder generar datos de ejemplo para probar la medida N50.
