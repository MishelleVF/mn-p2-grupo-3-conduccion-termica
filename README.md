# Análisis numérico de la conducción térmica unidimensional en barras metálicas

Repositorio del proyecto final del curso **CC2104 Métodos Numéricos**.  
El trabajo desarrolla el modelado numérico de la conducción térmica unidimensional en una barra metálica homogénea mediante métodos de diferencias finitas, con el propósito de analizar estabilidad, convergencia, precisión y robustez numérica.

## Descripción del proyecto

La conducción térmica unidimensional se modela mediante la ecuación del calor, considerando una barra metálica homogénea sometida a condiciones iniciales y de frontera definidas. El estudio compara el comportamiento de esquemas numéricos explícitos e implícitos para aproximar la evolución temporal de la temperatura.

Los métodos principales considerados son:

- Esquema explícito FTCS.
- Método implícito de Crank–Nicolson.
- Método iterativo de Gauss–Seidel para resolver los sistemas lineales asociados.
- Comparación con una solución analítica de referencia.
- Análisis de estabilidad mediante el parámetro de malla \(r\).

## Objetivo general

Modelar numéricamente la conducción térmica unidimensional en una barra metálica homogénea mediante métodos de diferencias finitas, con el propósito de analizar la estabilidad, convergencia y precisión de esquemas numéricos aplicados al estudio de la disipación de calor en componentes industriales alargados.

## Objetivos específicos

1. Formular e implementar el modelo numérico de la ecuación del calor unidimensional.
2. Evaluar el comportamiento del modelo bajo un conjunto base de parámetros espaciales, temporales y térmicos.
3. Analizar la robustez del modelo mediante la variación de parámetros relevantes, como la difusividad térmica, el tamaño de paso espacial y el tamaño de paso temporal.
4. Aplicar el modelo numérico al análisis de la disipación de calor en una barra metálica homogénea sometida a condiciones iniciales y de frontera definidas.

## Estructura del repositorio

```text
mn-p2-grupo-3-conduccion-termica/
│
├── main.mlx
├── ftcs_heat1d.m
├── crank_nicolson_heat1d.m
├── gauss_seidel_solver.m
├── exact_solution_heat1d.m
├── compute_errors.m
├── run_parameter_sweep.m
│
├── figures/
│   ├── temperatura_ftcs.png
│   ├── temperatura_crank_nicolson.png
│   ├── comparacion_metodos.png
│   ├── error_vs_dx.png
│   ├── error_vs_dt.png
│   └── estabilidad_r.png
│
├── tables/
│   ├── tabla_parametros.csv
│   ├── tabla_error_metodos.csv
│   └── tabla_costo_computacional.csv
│
├── paper/
│   └── MN_P2_Grupo_3.pdf
│
├── README.md
├── .gitignore
└── .gitattributes
```


## Integrantes:
* Marco Antonio Guerrero Ccompi
* Javier Ignacio Leon Olivares Sarmiento
* Marcelo Mateo Rosillo Rodriguez
* Enzo Marcelo Villanueva López
* Mishelle Stephany Villarreal Falcón

## Archivos principales
* main.mlx: archivo principal del proyecto en MATLAB Live Script. Ejecuta los experimentos, genera figuras, exporta tablas y organiza los resultados.
* ftcs_heat1d.m: implementación del esquema explícito FTCS.
* crank_nicolson_heat1d.m: implementación del método de Crank–Nicolson.
* gauss_seidel_solver.m: implementación del método iterativo de Gauss–Seidel.
* exact_solution_heat1d.m: cálculo de la solución analítica de referencia.
* compute_errors.m: cálculo de métricas de error entre soluciones numéricas y solución de referencia.
* run_parameter_sweep.m: ejecución de experimentos variando parámetros numéricos.

## Metodología computacional
- El flujo computacional seguido en el proyecto es:
1. Definición de parámetros físicos y numéricos del problema.
2. Discretización del dominio espacial y temporal.
3. Implementación del esquema FTCS.
4. Implementación del método de Crank–Nicolson.
5. Resolución del sistema lineal mediante Gauss–Seidel.
6. Comparación de resultados contra una solución analítica.
7. Cálculo de errores, estabilidad y costo computacional.
8. Generación de figuras y tablas para el paper final.

## Resultados esperados
- El proyecto busca obtener:
1.Evolución temporal de la temperatura en la barra.
2. Comparación entre FTCS y Crank–Nicolson.
3. Evaluación de estabilidad para distintos valores del parámetro r.
4. Análisis del error numérico frente a una solución exacta.
5. Comparación del costo computacional de los métodos.
6. Evidencia cuantitativa para la sección de Resultados y Discusión



