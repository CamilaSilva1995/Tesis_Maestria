Obtuvimos datos metagenomicos de Solena, estos son datos de
microrganismos en cultivos de fresa, para los cuales, inicialmente tenemos dos poblaciones, estas muestras estan etiquetadas como
como sanos y enfermos

## Scripts
<table class="default">
  <tr>
    <th scope="row">Nombre</th>
    <th>Explicación</th>
  </tr>
  <tr>
    <td><a href="https://github.com/CamilaSilva1995/Tesis_Maestria/blob/main/Analisis_Comparativo/Fresa_Solena/20230130_PreprocesamientoDatos.R">Preprocesamiento de Datos</a></td>
    <td>En este Scrip empezamos con un preprocesamiento de los datos, muestra como se descargaron los datos y como fue creado el archivo BIOM para ser poder leerlos desde R. Una vez se obtiene el archivo BIOM, este se varga en R como un objeto phyloseq. Y con este ya se procede a hacer una reconocimiento de los datos y primera observacion de medidas ecologicas como las diversidades alfa y beta. </td>
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
    <td><a href="">Rsss</a></td>
    <td>Esss</td>
  </tr>
  
</table>


## Reportes
Luego los scripts mostrados en la tabla anterior se convierten en los siguientes reportes, en los cuales se ordenan y se explican un poco mas los procesos, y sus respectivos resultados.

