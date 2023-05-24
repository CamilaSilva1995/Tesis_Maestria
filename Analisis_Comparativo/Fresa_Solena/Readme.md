# Análisis exploratorio de datos metagenómicos

El objetivo general de este análisis, es encontrar características diferenciadoras (marcadores funcionales) entre microbiomas de plantas sanas vs enfermas, usando datos metagenómicos. Para esto obtuvimos datos metagenómicos de <a href="https://solena.ag/home/us" >Solena</a>, estos son datos de microrganismos en cultivos de fresa, para los cuales, inicialmente tenemos dos poblaciones, estas muestras están etiquetadas como sanqs y enfermas.

Primero, se realizó un preprocesamiento de datos, tanto en Bash como en R; luego de esto se comenzó con el analisis exploratorio de los datos, el cual está principalmente dividido en tres partes: en primer lugar se realizó una exploración con diversidades alfa y beta, lo que llevo a un análisis estadístico con pruebas de hipótesis, en segundo lugar una visualización de correlación con redes y por ultimo una clasificación con machine learning.

A continuación, se muestra un resumen del análisis realizado, y a las concluciones obtenidas, explicando brevemente el contenido de cada script y markdown que contiene esta carpeta; empezando por los scripts:

## Scripts
<table class="default">
  <tr>
    <th scope="row">Nombre</th>
    <th>Explicación</th>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230130_PreprocesamientoDatos.R">Preprocesamiento de Datos</a></td>
    <td>En este Script empezamos con un preprocesamiento de los datos, muestra como se descargaron los datos y como fue creado el archivo BIOM para ser poder leerlos desde R. Una vez se obtiene el archivo BIOM, este se varga en R como un objeto phyloseq. Y con este ya se procede a hacer una reconocimiento de los datos y primera observacion de medidas ecologicas como las diversidades alfa y beta. </td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230213_DiversidadesAlfa%26Beta.R">Diversidades Alfa y Beta</a></td>
    <td>Aqui se puede ver la grafica de barras de abundancias de los datos y de los porcentajes de los datos en total, y tambein se empezo con distintas aglomeraciones de datos dependiendo del nivel taxonomico. Tambien se empiezan a ver las alfa y beta diversidad a los diferentes niveles taxonomicos.</td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230220_Funciones.R">Funciones</a></td>
    <td>Como mejora al anterior script, aqui se crean tres funciones automatizando el proceso de aglomeracion por niveles taxonomicos y sus respectivas graficas de barras de abundancia, alfa diversidad y beta diversidad.  </td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230227_Funciones%26Graficas.R">Funciones automatizadas</a></td>
    <td>Mejora de las funciones del script anterior, y creacion de todas las graficas, a niveles de Phylum, Genus, Specie y Familia, ceparado por Bacteria y Eukarya</td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230306_Potencia%26PruebaHipotesis.R">Potencia</a></td>
    <td>Analisis de potencia</td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230314_Redes.R">Redes</a></td>
    <td>En este Scrip se realizo un breve acercamiento a las redes simples y redes de coocurrencia, las redes siples usadas para la visualizacion, y las redes de coocurrencia tomando Fusarium como genero de interes en correlación. </td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230320_NuevosDatos.R">Nuevos datos</a></td>
    <td>Se obtuvo nuevos datos, los cuales estan divididos en tres nuevas categorias respecto a lugar donde se tomaron las muestras.</td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230321_NuevosDatosAll.R">Nuevos datos</a></td>
    <td>Tomando los nuevos datos del anterior script y los anteriores, se llego obtuvieron 5 categorias en total, tomando tanto el tipo de cultivo y si son muestras sanas o enfermas.</td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230327_PruebasOrden.R">PruebasOrden</a></td>
    <td>En este Scrip se realizo un reordenamiendo de las 5 categorias para poder visualizar mejor las graficas de barras.</td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230403_Actinobacteria.R">Actinobacteria</a></td>
    <td>Aqui vemos las graficas de barras de abundania y beta diversidad para Actinobacteria como filo de interes.</td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230403_Oomycota%26Fusarium.R">Fusarium</a></td>
    <td>Aqui vemos las graficas de barras de abundania y beta diversidad para Fusarium y Oomycota, como genero y filo de interes, respectivamente. </td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230410_PruebasdeHipotesis.R">Pruebas e Hipotesis</a></td>
    <td>Primer exploración con pruebas de hipotesis</td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230411_Preprocesamiento%26Normalizaci%C3%B3n.R">Preprocesamiento y Normalización</a></td>
    <td>Se realizó una normalización usando el paquete de R, "edgeR"</td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230419_Rarefaccion.R">Rarefaccion</a></td>
    <td>Se realizó una rarefacción de saturación, la cual nos indica que tenemos la cantidad necesaria de muestras para continuer con el analisis de los datos.</td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230425_Rarefaccion(Normalizacion).R">Rarefaccion</a></td>
    <td>Se realizo una normalización tipo rarefacción.</td>
  </tr>
   <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230427_PruebasdeHipotesisVarianzas.R">Prueba de hipótesis sobre varianzas</a></td>
     <td>Esss</td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230428_ML.ipynb">Clasificación con machine learning</a></td>
    <td>Esss</td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230502_NuevosDatos.R">Exproración de nuevos datos</a></td>
    <td>Esss</td>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230510_PruebaWilcoxon.R">Prueba de wilcoxon</a></td>
    <td>Esss</td>
  </tr>
</table>


## Reportes (R-MarkDown)
Luego los scripts mostrados en la tabla anterior se convierten en los siguientes reportes, en los cuales se ordenan y se explican un poco mas los procesos, y sus respectivos resultados.

<table class="default">
  <tr>
    <th scope="row">Nombre</th>
    <th>Explicación</th>
  </tr>
  <tr>
    <td><a href="https">Exploración</a></td>
    <td>En este documento podemos ver el preprocesamiento de los datos, una breve explicación del tipo de datos y tipo de formato en el cual se trabajaron, junto con una exploración inicial de los indices de diversidad alfa y beta, con diferentes medidas, a la totalidad de los datos; los cuales no mostraron una buena separacion de los datos entre muestras sanas y enfermas. Tambien se tiene una pequeña vizualización de correlacion con redes simples.</td>
  </tr>
  <tr>
    <td><a href="https">Exploración en conjuntos pequeños</a></td>
    <td>No se obtuvo una separacion visible de los datos en general, por lo tanto, en este reporte se realizaron subconjuntos por diferentes niveles taxonómicos, a los cuales tambien se visualizaron las diversidades alfa y beta.</td>
  </tr>
  <tr>
    <td><a href="https">Funciones automatizadas</a></td>
    <td>Se crearon funciones automatizadas para el análisis en conjuntos pequeños que se hizo en el reporte anterior.</td>
  </tr>
  <tr>
    <td><a href="https">Pruebas de Hipótesis</a></td>
    <td>Se realizaron distintas pruebas estadísticas, tanto al conjunto total de datos, sobre diferentes indices de diversidad; como a subconjuntos a distintos niveles taxonomicos de interes. Se realizo: prueba de Wilcoxon, pueba de hipótesis sobre las medias del indice Shannon y prueba de hipótesis sobre las varianzas de las medidas Chao1.</td>
  </tr>
  <tr>
    <td><a href="https">Aglomerados para Fusarium y Actinobacteria</a></td>
    <td>An...ia</td>
  </tr>
  <tr>
    <td><a href="https">Diversidades aglomerado 10%</a></td>
    <td>En...ción. </td>
  </tr>
  <tr>
    <td><a href="https">Redes de coocurrencia</a></td>
    <td>Se utilizaron los software (**AlnitaK** y MicNet) para crear las matrices de coocurrencia, tomando Fusarium como taxon principal,  y tambien de los datos en total.</td>
  </tr>
  <tr>
    <td><a href="https">Exploración para el segundo conjunto de datos (tres categorias) </a></td>
    <td>Tomando como referencia la primeta vis.</td>
  </tr>
  <tr>
    <td><a href="https">Exploración juntando los dos conjuntos de datos (cinco categorias) </a></td>
    <td>Eras.</td>
  </tr>
  <tr>
    <td><a href="https">Normalización </a></td>
    <td>Explorando diferntes tipos de rarefacción.</td>
  </tr>
  <tr>
    <td><a href="https">Exploración de datos normalizados</a></td>
    <td>Aqente. </td>
  </tr>
  <tr>
    <td><a href="https">Pruebas e Hipotesis</a></td>
    <td>Se realizaron diferentes tipos de pruebas de </td>
  </tr>
  <tr>
    <td><a href="https">Preprocesamiento y Normalización</a></td>
    <td>SR</td>
  </tr>
  <tr>
    <td><a href="https">Rarefaccion</a></<td>
    <td>Ses</td>
  </tr>
  <tr>
    <td><a href="https">Funciones automatizadas para el analisis en conjuntos pequeños, con datos normalizados</a></td>
    <td>Sn</td>
  </tr>
</table>

