# Conducción térmica 1D — Diferencias finitas y Crank–Nicolson

Código de la segunda entrega del proyecto de CC2104 (Métodos Numéricos).
Resuelve numéricamente la ecuación del calor 1D

    du/dt = alpha2 * d2u/dx2,   u(0,t)=u(L,t)=0,   u(x,0)=sin(pi*x/L)

y compara los esquemas **FTCS** (explícito) y **Crank–Nicolson** (implícito,
resuelto con **Gauss–Seidel**) contra la solución exacta.

## Estructura

```
matlab/
├── solucion_exacta.m       % solución analítica de referencia
├── ftcs.m                  % esquema explícito FTCS
├── gauss_seidel.m          % solver iterativo A*x=b
├── crank_nicolson_gs.m     % CN implícito + Gauss-Seidel por paso
├── calcular_error.m        % métricas de error vs exacta
├── run_experimentos.m      % Exp 1-4: figuras y tablas
└── main_conduccion_termica.m   % cuaderno principal (copiar a .mlx)
results/
├── figures/                % se generan al ejecutar
└── tables/                 % se generan al ejecutar
```

## Cómo ejecutar

1. Abrir MATLAB en la carpeta del repositorio.
2. Ejecutar el cuaderno principal:

   ```matlab
   run('matlab/main_conduccion_termica.m')
   ```

   o, para solo la batería de experimentos:

   ```matlab
   run('matlab/run_experimentos.m')
   ```

## Salidas generadas

Figuras (`results/figures/`):
- `fig1_perfiles_tiempos.png` — perfiles en varios tiempos (exacta/FTCS/CN).
- `fig2_tiempo_final.png` — comparación en el tiempo final.
- `fig3_estabilidad_ftcs.png` — FTCS estable (r≤0.5) vs inestable (r>0.5).
- `fig4_convergencia_Nx.png` — error L2 global vs Nx (log-log).

Tablas (`results/tables/`):
- `tabla1_comparacion_base.csv` — métricas por método (caso base).
- `tabla2_sensibilidad_Nx.csv` — refinamiento espacial.
- `tabla3_sensibilidad_r.csv` — variación del paso temporal.
- `tabla4_sensibilidad_alpha2.csv` — variación de la difusividad.

## Parámetros base

L=1, T=0.5, alpha2=0.01, Nx=50, r=0.25 (dt = r·dx²/alpha2).
