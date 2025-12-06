# Repositorio de código correspondiente al Trabajo Final de la Especialización en Explotación de Datos y Descubrimiento de Conocimiento de la U.B.A.

## Base de datos

Los archivos con el prefijo `nla22` o `nla2022` constituyen la base de datos (archivos con extensión `.csv.gz`) y su descripción (archivos con extensión `.txt`). Se trata de los resultados de los muestreos de la **National Lakes Assessmente 2022** de la E.P.A. de Estados Unidos, sobre los que este trabajo se construye.

## Código de otras fuentes
El archivo `cleanplot.pca.R` corresponde al siguiente libro:
Legendre, P. & L. Legendre. 2012. *Numerical ecology, 3rd English edition*

## Código

### Análisis

El archivo `main_script.r` es el que debe correrse si se desea realizar los análisis de nuevo (y generar los gráficos específicos de cada análisis). Éste incorpora y llama a `perform_analysis_new.R`. La totalidad de los análisis es de lenta ejecución.

En forma independiente se encuentra el script para realizar el análisis de outliers: `outlier_analysis.r`

### Resultados

Para generar las tablas de resumen debe correrse `results_script.r`, que tiene como input el archivo `all_results.RData`, que contiene el resultado final de la ejecución del archivo `main_script.r`. Es decir, que es posible obtener las tablas de resumen sin ejecutar de nuevo todos los análisis si así se lo desea.

Para generar los gráficos de resumen, comparando modelos, se debe ejecutar `plots_script.r`. Previamente, se tiene que haber corrido `results_script.r`, porque las tablas que genera son la entrada de este último script, y son necesarias para la generación de los gráficos.

### Licencia 

El código se encuentra disponible públicamente, y bajo la licencia GNU General Public License v3 (GPLv3).
