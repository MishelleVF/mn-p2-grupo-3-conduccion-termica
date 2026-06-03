function [x, t, U, r, iter_hist, gs_err_hist, info] = ...
         crank_nicolson_gs(L, T, alpha2, Nx, dt, f, tol, maxIter)
% CRANK_NICOLSON_GS  Esquema implicito de Crank-Nicolson para el calor 1D,
%   resolviendo el sistema lineal de cada paso temporal con Gauss-Seidel.
%
%   Proposito:
%       Resolver du/dt = alpha2 * d2u/dx2 con Dirichlet homogeneo mediante
%       Crank-Nicolson, que promedia las aproximaciones explicita e
%       implicita. Para los nodos interiores:
%           -(r/2) U_{i-1}^{n+1} + (1+r) U_i^{n+1} - (r/2) U_{i+1}^{n+1}
%         =  (r/2) U_{i-1}^{n}   + (1-r) U_i^{n}   + (r/2) U_{i+1}^{n}.
%       Esto equivale a A*u^{n+1} = B*u^{n}, resuelto en cada paso con
%       gauss_seidel.m.
%
%   Entradas:
%       L, T, alpha2, Nx, dt, f  - como en ftcs.m.
%       tol      - tolerancia para Gauss-Seidel.
%       maxIter  - iteraciones maximas de Gauss-Seidel por paso temporal.
%
%   Salidas:
%       x           - nodos espaciales (1 x (Nx+1)).
%       t           - instantes temporales (1 x (Nt+1)).
%       U           - matriz (Nt+1) x (Nx+1) de temperatura aproximada.
%       r           - parametro de malla r = alpha2*dt/dx^2.
%       iter_hist   - iteraciones de Gauss-Seidel en cada paso (Nt x 1).
%       gs_err_hist - error final de Gauss-Seidel en cada paso (Nt x 1).
%       info        - struct con metadata, incluye iter_prom (promedio GS).
%
%   Justificacion numerica:
%       A diferencia de FTCS, Crank-Nicolson es de segundo orden en tiempo
%       e INCONDICIONALMENTE estable (|G| <= 1 para todo r en la ecuacion
%       del calor lineal), pero exige resolver un sistema lineal por paso
%       porque los valores en n+1 estan acoplados. Comparar ambos esquemas
%       contra la solucion exacta permite contrastar precision, estabilidad
%       y costo (numero de iteraciones de Gauss-Seidel).

    % --- Validacion de entradas ---
    if L <= 0 || T <= 0 || dt <= 0 || Nx < 3
        error('cn_gs:params', 'Requiere L,T,dt>0 y Nx>=3.');
    end
    if ~isa(f, 'function_handle')
        error('cn_gs:f', 'f debe ser un handle de funcion.');
    end

    % --- Malla espacio-temporal ---
    dx = L / Nx;
    x  = linspace(0, L, Nx + 1);
    Nt = round(T / dt);
    t  = (0:Nt) * dt;
    r  = alpha2 * dt / dx^2;
    m  = Nx - 1;                      % numero de nodos interiores

    % --- Matrices tridiagonales A (implicita) y B (explicita) ---
    % Solo nodos interiores: con Dirichlet = 0 no hay aporte de frontera.
    A = diag((1 + r) * ones(m, 1)) ...
      + diag((-r/2)  * ones(m-1, 1),  1) ...
      + diag((-r/2)  * ones(m-1, 1), -1);
    B = diag((1 - r) * ones(m, 1)) ...
      + diag(( r/2)  * ones(m-1, 1),  1) ...
      + diag(( r/2)  * ones(m-1, 1), -1);

    % --- Condicion inicial y de frontera ---
    U = zeros(Nt + 1, Nx + 1);
    U(1, :)   = f(x);
    U(1, 1)   = 0;
    U(1, end) = 0;

    iter_hist   = zeros(Nt, 1);
    gs_err_hist = zeros(Nt, 1);

    % --- Avance temporal implicito ---
    for n = 1:Nt
        u_int = U(n, 2:Nx).';                 % interior actual (m x 1)
        rhs   = B * u_int;                    % lado derecho
        x0    = u_int;                        % buen punto de partida para GS
        [sol, it, e, conv] = gauss_seidel(A, rhs, x0, tol, maxIter);
        if ~conv
            warning('cn_gs:GS', ...
                'Gauss-Seidel no convergio en el paso %d (err=%.2e).', n, e);
        end
        U(n+1, 2:Nx) = sol.';
        U(n+1, 1)    = 0;                     % se conservan fronteras nulas
        U(n+1, end)  = 0;
        iter_hist(n)   = it;
        gs_err_hist(n) = e;
    end

    % --- Metadata ---
    info = struct('metodo', 'Crank-Nicolson + Gauss-Seidel', ...
                  'dx', dx, 'dt', dt, 'Nx', Nx, 'Nt', Nt, 'r', r, ...
                  'iter_prom', mean(iter_hist), 'tol', tol);
end
