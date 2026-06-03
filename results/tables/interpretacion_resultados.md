# Interpretación de resultados para el paper

Texto de referencia para la sección de Resultados y Discusión.  
Todos los valores provienen de la ejecución numérica; no se extrapolaron conclusiones no respaldadas por los datos.

---

## Figura 1 — Evolución temporal del perfil de temperatura

La Figura 1 muestra los perfiles de temperatura en cuatro instantes equidistribuidos en el intervalo `[0, T]`. La línea continua corresponde a la solución analítica, los círculos al esquema FTCS y los cuadrados a Crank–Nicolson. En todos los instantes, ambas soluciones numéricas se superponen visualmente con la solución exacta, lo que indica que el error de discretización es reducido para los parámetros base (`Nx = 50`, `r = 0.25`). El decaimiento progresivo de la amplitud refleja el comportamiento difusivo descrito por la ecuación del calor, con una tasa gobernada por la difusividad térmica `alpha2 = 0.01`.

---

## Figura 2 — Comparación en el tiempo final

El panel superior de la Figura 2 confronta los perfiles numéricos con la solución exacta en `t = T`. La coincidencia visual es casi perfecta, lo que puede inducir a concluir que ambos métodos tienen precisión equivalente. El panel inferior resuelve esta ambigüedad al representar el error absoluto puntual `|u_num - u_ex|`. Se observa que el error alcanza su máximo en torno a `x = 0.5`, posición en la que la condición inicial presenta su mayor amplitud, y que ambos métodos exhiben un perfil de error similar en forma aunque con magnitudes distintas. El análisis cuantitativo correspondiente se presenta en la Tabla Resumen.

---

## Figura 3 — Estabilidad de FTCS

La Figura 3 ilustra el criterio de estabilidad de Von Neumann para el esquema FTCS. El panel izquierdo corresponde al caso estable (`r = 0.40 ≤ 0.5`): el perfil numérico en `t = T` concuerda con la solución exacta y no se observan oscilaciones espurias. El panel derecho muestra el crecimiento de `max_x |u(x,t)|` en escala logarítmica para el caso inestable (`r = 0.75 > 0.5`, `Nx = 100`). La curva exhibe un crecimiento aproximadamente exponencial a partir de un umbral temporal, consistente con la amplificación de los modos de alta frecuencia predicha por el análisis espectral. Este comportamiento confirma que la condición `r ≤ 0.5` es necesaria para la estabilidad de FTCS y que su violación produce resultados sin validez física.

---

## Figura 4 — Convergencia espacial

La Figura 4 presenta el error L2 global en función del paso espacial `dx` para `Nx ∈ {25, 50, 100, 200}`, manteniendo `r = 0.25` fijo. La escala log-log permite visualizar directamente la pendiente de convergencia. Ambos métodos exhiben una pendiente aproximada de `p ≈ 2`, consistente con la discretización espacial centrada de segundo orden. La línea de referencia `p = 2` confirma el comportamiento asintótico. Los órdenes estimados de convergencia se reportan en la Tabla 2. Es relevante notar que la diferencia absoluta entre los errores de FTCS y Crank–Nicolson disminuye al refinar la malla, lo que sugiere que la ventaja relativa en precisión de un método sobre otro depende también de la resolución espacial.

---

## Tabla Resumen — Métricas del caso base

La Tabla Resumen consolida las métricas de error y el costo computacional de ambos métodos bajo los parámetros base (`Nx = 50`, `r = 0.25`, `alpha2 = 0.01`). Para este caso, FTCS presenta un error absoluto máximo y un RMSE ligeramente inferiores a los de Crank–Nicolson, lo que indica que bajo condiciones de estabilidad y parámetros base, la diferencia de precisión entre ambos esquemas no es sustancial. El tiempo de ejecución de FTCS es significativamente menor al de Crank–Nicolson, dado que el esquema explícito no requiere resolver un sistema lineal por paso temporal. Crank–Nicolson, en cambio, demanda iteraciones de Gauss–Seidel en cada paso, lo que eleva su costo computacional aunque garantiza estabilidad incondicional para cualquier valor de `r`.

---

## Tabla 2 — Sensibilidad al refinamiento espacial y orden de convergencia

La Tabla 2 reporta las métricas de error para `Nx ∈ {25, 50, 100, 200}` con `r = 0.25` constante. Los órdenes de convergencia estimados para ambos métodos convergen hacia `p ≈ 2` al refinar la malla, confirmando la convergencia de segundo orden en espacio. El primer valor de `orden` es `NaN` por ausencia de referencia anterior. Los errores relativos (`err_rel_L2_global`) son comparables entre FTCS y Crank–Nicolson en cada nivel de refinamiento, lo que sugiere que la superioridad de uno sobre otro en términos de precisión no es sistemática bajo estas condiciones. El costo computacional de Crank–Nicolson escala de forma más pronunciada con `Nx` debido al aumento del tamaño del sistema lineal.

---

## Tabla 3 — Sensibilidad al parámetro r (paso temporal)

La Tabla 3 evalúa el comportamiento de ambos métodos para `r ∈ {0.10, 0.25, 0.40, 0.50, 0.75}`. Para los casos estables (`r ≤ 0.5`), los errores de FTCS y Crank–Nicolson son consistentes y de magnitud similar. El error de FTCS presenta una ligera variación con `r`, mientras que el de Crank–Nicolson permanece prácticamente constante, lo que es coherente con su estabilidad incondicional. Para `r = 0.75`, FTCS viola el criterio de Von Neumann; en la columna `max_abs_u_FTCS` se registra el valor máximo de `|u|` al final de la simulación como indicador del crecimiento numérico. Los errores de FTCS para este caso se omiten (`NaN`) dado que no constituyen una medida válida de precisión. Crank–Nicolson mantiene errores reducidos incluso para `r = 0.75`, evidenciando su robustez ante pasos temporales grandes.

---

## Tabla 4 — Sensibilidad a la difusividad térmica (alpha2)

La Tabla 4 examina el efecto de la difusividad térmica `alpha2 ∈ {0.005, 0.01, 0.05, 0.10}` con `r = 0.25` y `Nx = 50` fijos. El paso temporal `dt` se ajusta proporcionalmente a `1/alpha2` para mantener `r` constante; por ello, mayor `alpha2` implica menor `dt` y mayor número de pasos temporales `Nt`. Los errores absolutos varían con `alpha2` porque la solución exacta cambia de escala (a mayor difusividad, el decaimiento es más rápido y la amplitud al tiempo final es menor), pero los errores relativos `err_rel_L2_global` resultan más comparables entre casos. El costo computacional de ambos métodos aumenta con `alpha2` debido al mayor número de pasos, siendo este aumento más pronunciado en Crank–Nicolson por la necesidad de resolver más sistemas lineales.
