# Análisis numérico de la conducción térmica unidimensional en barras metálicas

Repositorio del proyecto final del curso **CC2104 Métodos Numéricos**.  
El trabajo desarrolla el modelado numérico de la conducción térmica unidimensional en una barra metálica homogénea mediante esquemas de diferencias finitas, con el propósito de analizar estabilidad, convergencia, precisión y robustez numérica.

## Descripción del proyecto

La ecuación de calor unidimensional se resuelve numéricamente comparando un esquema explícito (FTCS) y uno implícito (Crank–Nicolson), ambos validados contra una solución analítica exacta. Se estudian la estabilidad condicional de FTCS (criterio `r ≤ 0.5`), la convergencia espacial de segundo orden y la sensibilidad a los parámetros de malla y a la difusividad térmica.

## Cómo ejecutar

Desde la raíz del repositorio, en MATLAB:

```matlab
addpath('matlab')
run('matlab/run_experimentos.m')
```

Esto genera automáticamente todas las figuras en `results/figures/` y todas las tablas en `results/tables/`.

Para la versión comentada paso a paso, ejecutar `matlab/main_conduccion_termica.m` (o la versión Live Script `.mlx` equivalente). Los archivos `.m` son la base reproducible del proyecto; los `.mlx` documentan la ejecución interactiva.

## Estructura del repositorio

```
mn-p2-grupo-3-conduccion-termica/
│
├── matlab/
│   ├── run_experimentos.m          ← script principal de experimentos
│   ├── main_conduccion_termica.m   ← cuaderno principal (caso base)
│   ├── ftcs.m                      ← esquema FTCS (explícito)
│   ├── crank_nicolson_gs.m         ← Crank-Nicolson + Gauss-Seidel (implícito)
│   ├── gauss_seidel.m              ← solucionador iterativo
│   ├── solucion_exacta.m           ← solución analítica de referencia
│   └── calcular_error.m            ← métricas de error (L∞, L2, RMSE, relativas)
│
├── results/
│   ├── figures/
│   │   ├── fig1_perfiles_tiempos.png      ← evolución temporal T(x,t)
│   │   ├── fig2_tiempo_final.png          ← perfiles finales + error absoluto
│   │   ├── fig3_estabilidad_ftcs.png      ← estabilidad FTCS (estable vs inestable)
│   │   └── fig4_convergencia_Nx.png       ← convergencia espacial L2 vs dx
│   └── tables/
│       ├── tabla_resumen_paper.csv        ← tabla principal para el paper
│       ├── tabla1_comparacion_base.csv    ← igual que tabla_resumen_paper
│       ├── tabla2_sensibilidad_Nx.csv     ← convergencia espacial + órdenes
│       ├── tabla3_sensibilidad_r.csv      ← sensibilidad al paso temporal r
│       ├── tabla4_sensibilidad_alpha2.csv ← sensibilidad a la difusividad
│       └── interpretacion_resultados.md  ← texto académico por figura/tabla
│
├── README.md
├── .gitignore
└── .gitattributes
```

## Figuras generadas

| Archivo | Contenido |
|---|---|
| `fig1_perfiles_tiempos.png` | Perfiles de temperatura en 4 instantes; exacta (línea), FTCS (círculos), CN (cuadrados) |
| `fig2_tiempo_final.png` | Subplot superior: perfiles en `t=T`; subplot inferior: error absoluto puntual |
| `fig3_estabilidad_ftcs.png` | Izquierda: FTCS estable `r=0.40`; derecha: crecimiento de `max|u|` para `r=0.75` (escala log) |
| `fig4_convergencia_Nx.png` | Error L2 global vs `dx` en escala log-log; referencia `p=2` |

## Tablas generadas

| Archivo | Contenido |
|---|---|
| `tabla_resumen_paper.csv` | Caso base: todas las métricas de error, tiempos, observaciones |
| `tabla2_sensibilidad_Nx.csv` | Refinamiento espacial: `Nx = [25, 50, 100, 200]`; incluye órdenes de convergencia |
| `tabla3_sensibilidad_r.csv` | Variación de `r = [0.10, 0.25, 0.40, 0.50, 0.75]`; caso `r=0.75` marcado como inestable |
| `tabla4_sensibilidad_alpha2.csv` | Variación de `alpha2` con `r` fijo; incluye tiempos de ejecución |

## Métricas de error implementadas (`calcular_error.m`)

- `max_abs` — error absoluto máximo (L∞)
- `L2_global` — norma L2 discreta espacio-tiempo
- `L2_final` — norma L2 discreta en el tiempo final
- `rel_max`, `rel_L2_global`, `rel_L2_final` — versiones relativas (normalizadas por la solución exacta)
- `RMSE_global`, `RMSE_final` — error cuadrático medio

## Integrantes

- Marco Antonio Guerrero Ccompi  
- Javier Ignacio Leon Olivares Sarmiento  
- Marcelo Mateo Rosillo Rodriguez  
- Enzo Marcelo Villanueva López  
- Mishelle Stephany Villarreal Falcón  

## Objetivos específicos

1. Formular e implementar el modelo numérico de la ecuación del calor unidimensional.  
2. Evaluar el comportamiento del modelo bajo un conjunto base de parámetros.  
3. Analizar la robustez del modelo variando la difusividad térmica, el paso espacial y el paso temporal.  
4. Aplicar el modelo al análisis de la disipación de calor en una barra metálica homogénea.
