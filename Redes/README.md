# Redes bayesianas (trabajo en curso)

Carpeta reservada para el modelado con redes bayesianas descrito en el capítulo 6 de la tesis.
Por ahora está vacía; el código se está desarrollando fuera del repositorio y se incorporará aquí.

## Objetivo

Modelar la estructura de dependencia probabilística entre el género fúngico *Fusarium*, el índice
de Shannon, los taxones más abundantes y el estado de salud de la planta, con un grafo acíclico
dirigido cuya estructura se aprende de los datos de `Data/fresa_solena/Data1/`.

## Plan

1. Construir la tabla de variables por muestra: abundancias relativas de los taxones más
   abundantes, abundancia de *Fusarium*, índice de Shannon y estado (`healthy` / `wilted`).
2. Aprender la estructura con el paquete `bnlearn` de R, por puntaje (BIC) o por restricciones.
3. Evaluar la robustez de los arcos con *bootstrap* y conservar los que aparecen en una proporción
   alta de las redes remuestreadas.
4. Medir la capacidad predictiva con validación cruzada de k particiones.

## Antecedentes en el repositorio

- Redes de coocurrencia con MicNet y Alnitak: `Analisis_Comparativo/Fresa_Solena/07_RedesCoocurrencia.Rmd`
  y `20230314_Redes.R`, con salidas en `Data/fresa_solena/Data1/Redes/`.
- Marco teórico: sección de redes bayesianas del capítulo 1 y sección 6.1 de la tesis.
