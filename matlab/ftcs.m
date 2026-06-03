function [x, t, U, r, info] = ftcs(L, T, alpha2, Nx, dt, f)
% FTCS  Esquema explicito (Forward Time, Centered Space) para el calor 1D.
%
%   Proposito:
%       Resolver du/dt = alpha2 * d2u/dx2 en [0,L] x [0,T] con condiciones
%       de Dirichlet homogeneas, mediante el esquema explicito
%           U_i^{n+1} = r*U_{i-1}^n + (1-2r)*U_i^n + r*U_{i+1}^n,
%       donde r = alpha2*dt/dx^2.
%
%   Entradas:
%       L       - longitud de la barra.
%       T       - tiempo final de simulacion.
%       alpha2  - difusividad termica (alpha^2).
%       Nx      - numero de subintervalos espaciales (=> Nx+1 nodos).
%       dt      - paso temporal. Para fijar un r objetivo o un Nt deseado:
%                   dt = r_obj*dx^2/alpha2   o   dt = T/Nt.
%       f       - handle de la condicion inicial, f(x). Ej: @(x) sin(pi*x/L).
%
%   Salidas:
%       x    - vector de nodos espaciales (1 x (Nx+1)).
%       t    - vector de instantes temporales (1 x (Nt+1)).
%       U    - matriz (Nt+1) x (Nx+1) con la temperatura aproximada.
%       r    - parametro de malla r = alpha2*dt/dx^2.
%       info - struct con metadata: metodo, dx, dt, Nx, Nt, r, estable.
%
%   Justificacion numerica:
%       FTCS aproxima du/dt por diferencia hacia adelante y d2u/dx2 por
%       diferencia centrada. El analisis de Von Neumann da el factor de
%       amplificacion G = 1 - 2r(1 - cos(theta)); el modo mas alto
%       (theta ~ pi) cumple |G| = |1 - 4r|, de modo que |G| <= 1 exige
%       0 <= r <= 1/2. Por eso el esquema es CONDICIONALMENTE estable:
%       si r > 1/2 los modos de alta frecuencia crecen y la solucion
%       diverge. La funcion NO se detiene si r > 1/2; devuelve el flag
%       info.estable para poder estudiar la inestabilidad explicitamente.

    % --- Validacion de entradas ---
    if L <= 0 || T <= 0 || dt <= 0 || Nx < 2
        error('ftcs:params', 'Parametros invalidos: requiere L,T,dt>0 y Nx>=2.');
    end
    if ~isa(f, 'function_handle')
        error('ftcs:f', 'f debe ser un handle, ej. @(x) sin(pi*x/L).');
    end

    % --- Malla espacio-temporal ---
    dx = L / Nx;
    x  = linspace(0, L, Nx + 1);     % nodos 1..Nx+1 (interiores: 2..Nx)
    Nt = round(T / dt);
    t  = (0:Nt) * dt;
    r  = alpha2 * dt / dx^2;
    estable = (r <= 0.5);            % condicion de estabilidad FTCS

    % --- Condicion inicial y de frontera ---
    U = zeros(Nt + 1, Nx + 1);
    U(1, :)   = f(x);
    U(1, 1)   = 0;                   % Dirichlet u(0,t)=0
    U(1, end) = 0;                   % Dirichlet u(L,t)=0

    % --- Avance temporal explicito (vectorizado en el espacio) ---
    in = 2:Nx;                       % indices de nodos interiores
    for n = 1:Nt
        U(n+1, in)  = r*U(n, in-1) + (1 - 2*r)*U(n, in) + r*U(n, in+1);
        U(n+1, 1)   = 0;             % se conservan las fronteras nulas
        U(n+1, end) = 0;
    end

    % --- Metadata del experimento ---
    info = struct('metodo', 'FTCS', 'dx', dx, 'dt', dt, ...
                  'Nx', Nx, 'Nt', Nt, 'r', r, 'estable', estable);
end
