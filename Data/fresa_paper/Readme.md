Estos datos estan publicos, y fueron tomados del articulo "(https://link.springer.com/article/10.1007/s00284-020-01948-x)", para llegar a una recreación del estudio que hacen en dicho articulo. por lo tanto se procedio a descargar los datos de NCBI, usando "SRA tools" con (fasta_dump), obteniendo asi los archivos .fasta de las muestras seleccionadas.

Luego para llegar a las tablas de abundancia se debe pasar por Qiime (dad y greengenes), ya que estos datos son de 16S.

Con referencia de el articulo y los datos encontrados en NCBI sobre las muestras usadas se creo el archivo de metadatos.



Para datos shotgun es preferible usar kraken, pero para estos datos probamos corriendo kraken (kraken2) para obtener los .report.

Con estas lineas corremos kraken para cada una de las muestras: (/files/kraken/db_kraken2 aqui la base de datos de kraken)

```{bash}
kraken2 --db /files/kraken/db_kraken2 -threads 14 --paired SRR8413967_1.fastq SRR8413967_2.fastq --output reports_fresa/SRR8413967.kraken --report reports_fresa/SRR8413967.report
```

Luego con estas lineas creamos el archivo .biom:

```{bash}
kraken-biom SRR8413967.report  SRR8413968.report  SRR8413969.report  SRR8413970.report --fmt json -o fresa_paper.biom
```
Todo esto activando el ambiente metagenomics previamente (conda activate metagenomics)

Luego de esto ya es posible usar los datos para el analisis exploratorio en R.
