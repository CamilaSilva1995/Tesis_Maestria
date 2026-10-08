# Taller práctico: pruebas de hipótesis con datos de microbioma

Taller en Python, pensado para Google Colab y dirigido a quien quiera entender las pruebas
estadísticas que se usan para comparar grupos de muestras metagenómicas. Está organizado en siete
episodios, al estilo de las lecciones de The Carpentries: cada uno abre con preguntas y objetivos,
sigue con la explicación paso a paso y cierra con un ejercicio y sus puntos clave. Dura unas 3 horas
con descanso (2 h de explicación y 56 min de ejercicios) y puede darse en dos sesiones.

Hay además un **taller rápido**, una versión resumida de 57 minutos en cuatro episodios, con los
mismos datos y los mismos resultados, para cuando solo se dispone de una hora.

Cada prueba se aplica a dos casos:

- **Datos reales**, los de la tesis: 53 metagenomas de rizósfera de fresa, 35 de plantas saludables
  y 18 de plantas no saludables.
- **Datos simulados** dentro del cuaderno: 60 muestras (30 por grupo) que cumplen los supuestos de
  las pruebas y en las que se conoce la verdad, porque la elige quien simula.

Temas: qué es una prueba de hipótesis y el valor p (con permutaciones), prueba t de Student y de
Welch con sus supuestos, Mann-Whitney y Wilcoxon de rangos con signo, proporción de un género por
grupo (*Fusarium* en los datos reales), y PERMANOVA con PERMDISP y una ordenación. Tiene siete
ejercicios, uno por episodio, con la solución escondida en un desplegable.

## Conocimientos previos

El taller empieza donde terminan estas lecciones de The Carpentries, que el cuaderno cita como punto
de partida recomendado:

- [Data Processing and Visualization for Metagenomics](https://carpentries-lab.github.io/metagenomics-analysis/)
  (Carpentries Lab): cómo se llega de las lecturas a la tabla de conteos y a las gráficas de diversidad.
- [Plotting and Programming in Python](https://swcarpentry.github.io/python-novice-gapminder/)
  (Software Carpentry): variables, funciones, bucles y gráficos.
- [Data Analysis and Visualization in Python for Ecologists](https://datacarpentry.github.io/python-ecology-lesson/)
  (Data Carpentry), con versión en español,
  [Análisis y visualización de datos usando Python](https://datacarpentry.github.io/python-ecology-lesson-es/):
  tablas de pandas.

## Abrir en Google Colab

1. Entra a [colab.research.google.com](https://colab.research.google.com) y elige
   **Archivo → Abrir notebook → GitHub**.
2. Pega la dirección del repositorio, `CamilaSilva1995/Tesis_Maestria`, y abre
   `Taller_Practico/Taller_Practico_Analisis_estadistico_de_Datos_Metagenomicos.ipynb`.
3. Ejecuta las celdas en orden. Los datos reales se descargan solos desde GitHub y los simulados se
   generan en el propio cuaderno; no hay que subir nada ni instalar paquetes, porque solo usa numpy,
   pandas, scipy y matplotlib.

También puede abrirse en Jupyter o VS Code desde esta carpeta; en ese caso los datos se leen de
`datos/`.

## Archivos

| Archivo | Contenido |
|---|---|
| `Taller_Practico_Analisis_estadistico_de_Datos_Metagenomicos.ipynb` | El cuaderno del taller, con las salidas ya ejecutadas para que los estudiantes vean el resultado esperado. |
| `construir_cuaderno.py` | Script que genera el cuaderno episodio por episodio y calcula el tiempo de cada uno. Para cambiar el texto, el código o los ejercicios, edita este archivo y vuelve a ejecutarlo. |
| `Taller_Rapido_Pruebas_de_hipotesis_con_Datos_Metagenomicos.ipynb` | El taller rápido: la versión resumida de una hora, también con las salidas ya ejecutadas. |
| `construir_taller_rapido.py` | Script que genera el taller rápido. Se detiene con un error si el contenido pasa de 60 minutos. |
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
# ejecutarlo para guardar las salidas
python3 -c "
import nbformat; from nbconvert.preprocessors import ExecutePreprocessor
nb = nbformat.read('Taller_Practico_Analisis_estadistico_de_Datos_Metagenomicos.ipynb', as_version=4)
ExecutePreprocessor(timeout=600).preprocess(nb, {'metadata': {'path': '.'}})
nbformat.write(nb, 'Taller_Practico_Analisis_estadistico_de_Datos_Metagenomicos.ipynb')"
```

Al terminar, `construir_cuaderno.py` imprime de dónde sale el tiempo de cada episodio. El taller
rápido se regenera igual, con `construir_taller_rapido.py` y el nombre de su cuaderno.

## Resultados que obtiene el taller

**Datos reales.** Coinciden con el capítulo 4 de la tesis: Welch p = 0.150, Mann-Whitney U = 376 y
p = 0.259, prueba F p < 0.001, PERMANOVA R² ≈ 0.02 y p ≈ 0.28 sobre géneros. Además, PERMDISP da
p ≈ 0.005, lo que confirma que las plantas no saludables son más heterogéneas también en composición.

**Datos simulados** (semilla 2026). Los supuestos se cumplen (Shapiro-Wilk p = 0.84 y 0.19, prueba F
p = 0.11) y las pruebas coinciden: permutación p = 0.006, Student p = 0.0046, Welch p = 0.0047,
Mann-Whitney p = 0.003. El género simulado como patógeno se detecta (U = 0), el PERMANOVA es
significativo (R² ≈ 0.28, p = 0.001) y PERMDISP no lo es (p ≈ 0.65). Los datos simulados dependen del
generador de números aleatorios de numpy: con la semilla fija se repiten, pero una versión futura de
numpy podría cambiar los últimos decimales.

## Guía para la sesión

| Inicio | Episodio | Explicación | Ejercicio |
|---|---|---|---|
| 0:00 | 1. Los datos: un caso real y un caso ideal | 25 min | 8 min |
| 0:33 | 2. La lógica de una prueba de hipótesis | 15 min | 6 min |
| 0:54 | 3. Comparar medias: la prueba t de Student y la de Welch | 20 min | 8 min |
| 1:22 | Descanso | 10 min | |
| 1:32 | 4. Pruebas basadas en rangos: Mann-Whitney y Wilcoxon | 15 min | 6 min |
| 1:53 | 5. Un género de interés: ¿es más abundante en un grupo? | 10 min | 10 min |
| 2:13 | 6. Comparar la composición completa: PERMANOVA y PERMDISP | 25 min | 10 min |
| 2:48 | 7. Cierre: reportar y comparar los dos casos | 10 min | 8 min |
| 3:06 | Fin | | |

**Cómo se calcula el tiempo.** La explicación de cada episodio se estima a partir de su contenido:
120 palabras de texto por minuto, 2 minutos por celda de código (leer los comentarios, ejecutarla y
revisar la salida), 1 minuto por figura y 1 minuto por cada pregunta «Para pensar»; el resultado se
redondea hacia arriba al múltiplo de 5 para dejar margen a las preguntas. Cada ejercicio trae su
propio tiempo. Las constantes están al inicio de `construir_cuaderno.py`: si cambias el contenido, el
índice del cuaderno se recalcula solo (esta tabla hay que actualizarla a mano).

**Cómo llevar cada episodio.** Leer las preguntas, explicar qué se hace y por qué, ejecutar primero
el caso simulado, hacer la pregunta «Para pensar» antes de ejecutar el caso real, dejar el ejercicio
en parejas y cerrar con los puntos clave.

**Si hay menos tiempo.** En dos sesiones: episodios 1 a 3 (1 h 22 min) y 4 a 7 (1 h 34 min). Para una
sola sesión de una hora, usar el taller rápido.

### Taller rápido (57 min)

| Inicio | Episodio | Explicación | Ejercicio |
|---|---|---|---|
| 0:00 | 1. Dos casos y una pregunta | 10 min | |
| 0:10 | 2. ¿Diferencia real o azar? El valor p | 10 min | 4 min |
| 0:24 | 3. Elegir la prueba: supuestos, Welch y Mann-Whitney | 10 min | 4 min |
| 0:38 | 4. Comparar comunidades completas: PERMANOVA y PERMDISP | 15 min | 4 min |
| 0:57 | Fin | | |

Usa la misma regla de tiempo y las mismas semillas que el taller completo, así que los números
coinciden. Deja fuera Wilcoxon de rangos con signo, la comparación de un género (*Fusarium*) con las
comparaciones múltiples, la construcción paso a paso del pseudo-F y la redacción del resultado; el
cuaderno remite al taller completo para esos temas.
