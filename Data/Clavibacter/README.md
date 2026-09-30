# Genomas de Clavibacter

Genomas de referencia del género *Clavibacter*, bacterias fitopatógenas de interés para Solena,
usados como material de prueba para el simulador de lecturas PyMetaSeem (`Generador_de_reads/`) y
para la métrica N50 (`N50/`).

| Subcarpeta | Contenido |
|---|---|
| `all-c-genomes/` | Genomas completos en formato FASTA descargados del NCBI. `all-c-genomes.zip` es el paquete original. |
| `all-c-genomes(cortados)/` | Los mismos genomas con el encabezado FASTA simplificado, que es el formato que espera PyMetaSeem. |
| `new-c-genomes(cortados)/` | Genomas adicionales incorporados después, ya recortados. |

El recorte del encabezado se hace con el script `cutout.py` de `Comparacion_taxonomica/` o
directamente con `sed`, como se indica al inicio del cuaderno `PyMetaSeem_Cami.ipynb`.
