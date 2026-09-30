# Comunidad simulada de diez genomas

Genomas de referencia de una comunidad microbiana simulada, usados como entrada del simulador de
lecturas PyMetaSeem (`Generador_de_reads/`). La composición imita la de los estándares comerciales
de comunidades simuladas: ocho bacterias y dos levaduras.

| Archivo | Organismo |
|---|---|
| `0. Listeria monocytogenes.csv` | *Listeria monocytogenes* |
| `1. Pseudomonas aeruginosa.csv` | *Pseudomonas aeruginosa* |
| `2. Bacillus subtilis.csv` | *Bacillus subtilis* |
| `3. Escherichia coli.csv` | *Escherichia coli* |
| `4. Salmonella enterica.csv` | *Salmonella enterica* |
| `5. Lactobacillus fermentum.csv` | *Lactobacillus fermentum* |
| `6. Enterococcus faecalis.csv` | *Enterococcus faecalis* |
| `7. Staphylococcus aureus.csv` | *Staphylococcus aureus* |
| `8. Saccharomyces cerevisiae.csv` | *Saccharomyces cerevisiae* |
| `9. Cryptococcus neoformans.csv` | *Cryptococcus neoformans* |

Los `.csv` son las tablas de contigs con sus longitudes por genoma.

| Subcarpeta | Contenido |
|---|---|
| `genomes/` | Genomas completos en FASTA. |
| `cut_head/` | Los mismos genomas con el encabezado FASTA simplificado, más `table.txt`, el perfil taxonómico con el número de lecturas que PyMetaSeem debe generar por genoma. |
| `tarea_Fonty/` | Ejercicio de clase que dio origen a esta comunidad. |

El cuaderno `Generador_de_reads/PyMetaSeem_Cf.ipynb` lee `cut_head/table.txt` y genera un archivo
FASTQ por genoma con la cantidad de lecturas indicada.
