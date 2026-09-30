# Análisis comparativo

Análisis exploratorio y estadístico de comunidades microbianas a partir de tablas de abundancia
generadas con Kraken. En todos los casos el flujo es el mismo: cargar el archivo BIOM como objeto
phyloseq en R, filtrar por calidad, calcular diversidades alfa y beta, explorar la composición por
nivel taxonómico y, cuando corresponde, contrastar diferencias con pruebas de hipótesis.

| Subcarpeta | Datos | Contenido |
|---|---|---|
| [`Fresa_Solena/`](Fresa_Solena/) | `Data/fresa_solena/` | Análisis principal de la tesis: rizósfera de fresa, plantas saludables frente a no saludables. Cuadernos R Markdown numerados, scripts de R y scripts que regeneran las figuras de la tesis. |
| [`Fresa_paper/`](Fresa_paper/) | `Data/fresa_paper/` | Réplica del análisis sobre los datos públicos de Yang et al. (2020), fresa con y sin oídio. |
| [`Clavi_Solena/`](Clavi_Solena/) | `Data/solena/` | Primer ejercicio con datos de Solena: *Clavibacter michiganensis* en cultivos de chile, maíz y tomate. |

El análisis de `Fresa_Solena/` es el que se reporta en los capítulos 3 a 5 de la tesis y en el
artículo de `Articulo_Popayan/`.
