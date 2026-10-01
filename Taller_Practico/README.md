# Taller práctico: pruebas de hipótesis con datos de microbioma

Clase práctica de 50 minutos, en Python y pensada para Google Colab, que enseña a comparar dos
grupos de muestras metagenómicas con los datos reales de la tesis: 53 metagenomas de rizósfera
de fresa, 35 de plantas saludables y 18 de plantas no saludables.

Temas: qué es una prueba de hipótesis y el valor p (con una simulación por permutaciones), prueba
t y prueba de Welch, Mann-Whitney y Wilcoxon de rangos con signo, proporción de *Fusarium* por
grupo, y PERMANOVA con PERMDISP. Tiene siete ejercicios con la solución escondida en un
desplegable.

## Abrir en Google Colab

1. Entra a [colab.research.google.com](https://colab.research.google.com) y elige
   **Archivo → Abrir notebook → GitHub**.
2. Pega la dirección del repositorio, `CamilaSilva1995/Tesis_Maestria`, y abre
   `Taller_Practico/Taller_Practico_Analisis_estadistico_de_Datos_Metagenomicos.ipynb`.
3. Ejecuta las celdas en orden. Los datos se descargan solos desde GitHub; no hay que subir nada
   ni instalar paquetes, porque solo usa numpy, pandas, scipy y matplotlib.

También puede abrirse en Jupyter o VS Code desde esta carpeta; en ese caso los datos se leen de
`datos/`.

## Archivos

| Archivo | Contenido |
|---|---|
| `Taller_Practico_Analisis_estadistico_de_Datos_Metagenomicos.ipynb` | El cuaderno de la clase, con las salidas ya ejecutadas para que los estudiantes vean el resultado esperado. |
| `construir_cuaderno.py` | Script que genera el cuaderno celda por celda. Para cambiar el texto o los ejercicios, edita este archivo y vuelve a ejecutarlo. |
| `datos/metadatos.csv` | Una fila por muestra (58): grupo y número de lecturas clasificadas. El cuaderno aplica el filtro de calidad que deja 53. |
| `datos/diversidad_alfa.csv` | Riqueza observada, Chao1, Shannon y Simpson por muestra, calculados con phyloseq sobre la tabla completa, igual que en la tesis. |
| `datos/conteos_por_genero.csv` | Lecturas por género y muestra: 1 795 géneros, con su reino y filo. |
| `datos/conteos_por_filo.csv` | Lecturas por filo y muestra, por si se quiere trabajar a ese nivel. |
| `original/Practica_microbioma_fresa.ipynb` | Primera versión del cuaderno, más extensa y con los datos incrustados. Se conserva como referencia. |

Los tres CSV se generaron a partir de `Data/fresa_solena/Data1/fresa_kraken.biom` y
`metadata.csv`, los mismos archivos que usa la tesis.

## Regenerar el cuaderno

```bash
cd Taller_Practico
python3 construir_cuaderno.py
# opcional: ejecutarlo para guardar las salidas
python3 -c "
import nbformat; from nbconvert.preprocessors import ExecutePreprocessor
nb = nbformat.read('Taller_Practico_Analisis_estadistico_de_Datos_Metagenomicos.ipynb', as_version=4)
ExecutePreprocessor(timeout=600).preprocess(nb, {'metadata': {'path': '.'}})
nbformat.write(nb, 'Taller_Practico_Analisis_estadistico_de_Datos_Metagenomicos.ipynb')"
```

## Resultados que obtiene la clase

Coinciden con el capítulo 4 de la tesis: Welch p = 0.150, Mann-Whitney U = 376 y p = 0.259,
prueba F p < 0.001, PERMANOVA R² ≈ 0.02 y p ≈ 0.28 sobre géneros. Además, PERMDISP da p ≈ 0.005,
lo que confirma que las plantas no saludables son más heterogéneas también en composición.

## Guía para la sesión

| Minutos | Sección | Qué hacer |
|---|---|---|
| 0–8 | Datos | Cargar, filtrar y unir por identificador. Ejercicio 1 en parejas. |
| 8–14 | Prueba de hipótesis | Correr la simulación y discutir la figura antes de leer el texto. Ejercicio 2 oral. |
| 14–22 | t y Welch | Mostrar cómo el supuesto de varianzas cambia el valor p. Ejercicio 3. |
| 22–29 | Mann-Whitney y Wilcoxon | Insistir en independientes frente a pareados. Ejercicio 4 oral. |
| 29–37 | *Fusarium* | Discutir el denominador y la lección del emparejamiento de tablas. Ejercicio 5. |
| 37–47 | PERMANOVA | Construirlo a mano paso a paso; PERMDISP como advertencia. Ejercicio 6 si hay tiempo. |
| 47–50 | Cierre | Tabla resumen. Ejercicio 7 como tarea. |
