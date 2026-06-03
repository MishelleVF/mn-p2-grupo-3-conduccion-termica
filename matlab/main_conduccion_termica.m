%% Modelado numerico de la conduccion termica unidimensional
%  mediante diferencias finitas y Crank-Nicolson
%
%  Cuaderno principal de la segunda entrega (CC2104 - Metodos Numericos).
%  Pensado para copiarse a main_conduccion_termica.mlx: cada bloque %%
%  corresponde a una seccion del Live Script.
%
%  Funciones del proyecto (carpeta matlab/):
%     ftcs.m, crank_nicolson_gs.m, gauss_seidel.m,
%     solucion_exacta.m, calcular_error.m, run_experimentos.m
%
%% 1. Problema y modelo
%  Se estudia la conduccion termica 1D en una barra metalica homogenea:
%
%     du/dt = alpha2 * d2u/dx2,   x en [0,L],  t en [0,T]
%
%  Condiciones de frontera (Dirichlet homogeneas):  u(0,t)=u(L,t)=0
%  Condicion inicial:                                u(x,0)=sin(pi*x/L)
%  Solucion exacta de referencia:
%     u(x,t) = exp(-alpha2*(pi/L)^2 * t) * sin(pi*x/L)
%
%  Se comparan dos esquemas de diferencias finitas:
%     - FTCS           : explicito, condicionalmente estable (r <= 1/2).
%     - Crank-Nicolson : implicito, incondicionalmente estable; su sistema
%                        lineal por paso se resuelve con Gauss-Seidel.

clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir); script_dir = pwd; end
addpath(script_dir);

%% 2. Parametros base
L      = 1;          % longitud de la barra
T      = 0.5;        % tiempo final
alpha2 = 0.01;       % difusividad termica
Nx     = 50;         % subintervalos espaciales (Nx+1 nodos)
f      = @(x) sin(pi * x / L);   % condicion inicial
tol    = 1e-10;
maxIter = 10000;

dx = L / Nx;
dt = 0.25 * dx^2 / alpha2;        % r = 0.25 -> FTCS estable
fprintf('dx=%.4f  dt=%.5f  r=%.3f\n', dx, dt, alpha2*dt/dx^2);

%% 3. Ejecucion de los metodos (OE1)
[xf, tf, Uf, rf, infoF] = ftcs(L, T, alpha2, Nx, dt, f);
[xc, tc, Uc, rc, itC, gsErrC, infoC] = ...
        crank_nicolson_gs(L, T, alpha2, Nx, dt, f, tol, maxIter);
Ue = solucion_exacta(xf, tf, L, alpha2);

%% 4. Evaluacion del comportamiento (OE2): errores
errF = calcular_error(Uf, Ue, xf, tf);
errC = calcular_error(Uc, Ue, xc, tc);

TablaErr = table( {'FTCS';'Crank-Nicolson'}, ...
                  [errF.max_abs;   errC.max_abs], ...
                  [errF.L2_global; errC.L2_global], ...
                  [errF.L2_final;  errC.L2_final], ...
                  [NaN; infoC.iter_prom], ...
    'VariableNames', {'Metodo','err_max','err_L2_global','err_L2_final','GS_iter_prom'});
disp(TablaErr);

%% 5. Graficos principales
% 5.1 Evolucion espacio-temporal (CN)
figure('Color','w');
surf(xc, tc, Uc, 'EdgeColor','none'); view(135,30);
xlabel('x'); ylabel('t'); zlabel('u(x,t)');
title('Evolucion de la temperatura (Crank-Nicolson)'); colorbar;

% 5.2 Comparacion en el tiempo final
figure('Color','w'); hold on; grid on;
plot(xf, Ue(end,:), 'k-',  'LineWidth', 2);
plot(xf, Uf(end,:), 'bo--');
plot(xc, Uc(end,:), 'rs--');
xlabel('x'); ylabel(sprintf('u(x, %.3f)', T));
title('Tiempo final: exacta vs FTCS vs Crank-Nicolson');
legend('Exacta','FTCS','Crank-Nicolson');

%% 6. Experimentos completos (OE2, OE3, OE4)
%  Estabilidad de FTCS y estudio de sensibilidad (Nx, r, alpha2).
%  Genera y guarda todas las figuras y tablas en results/.
run_experimentos;

%% 7. Nota de alineacion con los objetivos
%  OE1  -> sec.3: implementacion de FTCS, CN y solucion exacta.
%  OE2  -> sec.4-5 y Exp1-2: comportamiento y comparacion cuantitativa base.
%  OE3  -> Exp3 (estabilidad FTCS) y Exp4 (sensibilidad a Nx, r, alpha2).
%  OE4  -> perfiles de disipacion del calor en la barra (figuras de perfil).
