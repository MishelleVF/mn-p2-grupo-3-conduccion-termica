%% run_experimentos.m
% Ejecuta los experimentos de la segunda entrega del proyecto:
%   "Modelado numerico de la conduccion termica unidimensional mediante
%    diferencias finitas y Crank-Nicolson".
%
% Genera figuras (results/figures) y tablas (results/tables) que responden
% directamente a los objetivos especificos:
%   OE1 Formular e implementar el modelo numerico (FTCS, CN, exacta).
%   OE2 Evaluar el modelo bajo un conjunto base de parametros.
%   OE3 Analizar la robustez variando alpha2, dx (Nx) y dt (r).
%   OE4 Aplicar el modelo a la disipacion de calor en una barra alargada.
%
% Para ejecutar desde la raiz del repositorio:
%   addpath('matlab')
%   run('matlab/run_experimentos.m')
%
% NO inventa resultados: todo proviene de la ejecucion numerica.

clear; clc; close all;

%% -------- Configuracion de rutas --------
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
addpath(script_dir);
repo_root = fileparts(script_dir);
if isempty(repo_root)
    repo_root = pwd;
end
fig_dir = fullfile(repo_root, 'results', 'figures');
tab_dir = fullfile(repo_root, 'results', 'tables');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
if ~exist(tab_dir, 'dir'); mkdir(tab_dir); end

% Caracteres con tilde (char() garantiza portabilidad en cualquier MATLAB)
o_ac = char(243);   % o con tilde: o
e_ac = char(233);   % e con tilde: e
i_ac = char(237);   % i con tilde: i
a_ac = char(225);   % a con tilde: a

%% -------- Parametros base (normalizados) --------
L       = 1;          % longitud de la barra
T       = 0.5;        % tiempo final
alpha2  = 0.01;       % difusividad termica
Nx      = 50;         % subintervalos espaciales
f       = @(x) sin(pi * x / L);   % condicion inicial compatible con la exacta
tol     = 1e-10;      % tolerancia de Gauss-Seidel
maxIter = 10000;      % iteraciones maximas de Gauss-Seidel

dx     = L / Nx;
r_base = 0.25;                       % r objetivo estable para FTCS
dt     = r_base * dx^2 / alpha2;     % => r = 0.25

fprintf('Parametros base: L=%.2f T=%.2f alpha2=%.4f Nx=%d dx=%.4f dt=%.5f r=%.3f\n', ...
        L, T, alpha2, Nx, dx, dt, alpha2*dt/dx^2);

%% ===================================================================
%  EXPERIMENTO 1: CASO BASE  (OE1, OE2, OE4)
%  Ejecuta FTCS y Crank-Nicolson, compara con la solucion exacta.
%  ===================================================================
[xf, tf, Uf, rf, infoF] = ftcs(L, T, alpha2, Nx, dt, f);
[xc, tc, Uc, rc, itC, gsErrC, infoC] = ...
        crank_nicolson_gs(L, T, alpha2, Nx, dt, f, tol, maxIter);
Ue = solucion_exacta(xf, tf, L, alpha2);

errF = calcular_error(Uf, Ue, xf, tf);
errC = calcular_error(Uc, Ue, xc, tc);

% ---- Figura 1: perfiles de temperatura en distintos tiempos ----
% Leyenda en dos grupos: metodo (trazos negros dummy) + instante (colores)
idx = unique(round(linspace(1, numel(tf), 4)));
fig1 = figure('Color', 'w', 'Position', [100 100 820 520]);
hold on; grid on;
colores = lines(numel(idx));
h_time = gobjects(numel(idx), 1);
for k = 1:numel(idx)
    n = idx(k);
    h_time(k) = plot(xf, Ue(n,:), '-', 'Color', colores(k,:), 'LineWidth', 1.8);
    plot(xf, Uf(n,:), 'o', 'Color', colores(k,:), 'MarkerSize', 4, 'LineStyle', 'none');
    plot(xc, Uc(n,:), 's', 'Color', colores(k,:), 'MarkerSize', 4, 'LineStyle', 'none');
end
hE = plot(NaN, NaN, 'k-',  'LineWidth', 1.8);
hF = plot(NaN, NaN, 'ko',  'MarkerSize', 5, 'LineStyle', 'none');
hC = plot(NaN, NaN, 'ks',  'MarkerSize', 5, 'LineStyle', 'none');

xlabel(['Posici' o_ac 'n x']); ylabel('Temperatura u(x,t)');
title(['Evoluci' o_ac 'n temporal del perfil de temperatura en la barra'], ...
      'FontWeight', 'normal');
t_labels = arrayfun(@(n) sprintf('t = %.3f', tf(n)), idx, 'UniformOutput', false);
legend([hE; hF; hC; h_time], ...
       [{'Exacta', 'FTCS', 'Crank-Nicolson'}, t_labels], ...
       'Location', 'northeast', 'FontSize', 8);
saveas(fig1, fullfile(fig_dir, 'fig1_perfiles_tiempos.png'));

% ---- Figura 2: comparacion en el tiempo final + error absoluto puntual ----
% Subplot superior: perfiles exacto vs FTCS vs CN en t = T.
% Subplot inferior: error absoluto puntual |u_num - u_exacta| en t = T.
fig2 = figure('Color', 'w', 'Position', [100 100 720 560]);

subplot(2,1,1);
hold on; grid on;
plot(xf, Ue(end,:), 'k-',   'LineWidth', 2,  'DisplayName', 'Exacta');
plot(xf, Uf(end,:), 'bo--', 'MarkerSize', 4, 'DisplayName', 'FTCS');
plot(xc, Uc(end,:), 'rs--', 'MarkerSize', 4, 'DisplayName', 'Crank-Nicolson');
ylabel(sprintf('u(x, t = %.3f)', T));
title(['Comparaci' o_ac 'n de perfiles en el tiempo final'], 'FontWeight', 'normal');
legend('Location', 'northeast', 'FontSize', 9);

subplot(2,1,2);
hold on; grid on;
err_abs_F = abs(Uf(end,:) - Ue(end,:));
err_abs_C = abs(Uc(end,:) - Ue(end,:));
plot(xf, err_abs_F, 'b-', 'LineWidth', 1.4, 'DisplayName', 'FTCS');
plot(xc, err_abs_C, 'r-', 'LineWidth', 1.4, 'DisplayName', 'Crank-Nicolson');
xlabel(['Posici' o_ac 'n x']); ylabel('|u_{num} - u_{ex}|');
title('Error absoluto puntual en el tiempo final', 'FontWeight', 'normal');
legend('Location', 'north', 'FontSize', 9);

saveas(fig2, fullfile(fig_dir, 'fig2_tiempo_final.png'));

fprintf('\n[Exp1] FTCS: err_max=%.3e  L2_final=%.3e  rel_L2=%.3e  RMSE=%.3e\n', ...
        errF.max_abs, errF.L2_final, errF.rel_L2_global, errF.RMSE_global);
fprintf('[Exp1] CN  : err_max=%.3e  L2_final=%.3e  rel_L2=%.3e  RMSE=%.3e  GS_iter=%.1f\n', ...
        errC.max_abs, errC.L2_final, errC.rel_L2_global, errC.RMSE_global, infoC.iter_prom);

%% ===================================================================
%  EXPERIMENTO 2: TABLA COMPARATIVA Y TABLA RESUMEN  (OE2)
%  Tabla principal con todas las metricas de error para el paper.
%  ===================================================================
rep = 5;
tic; for k=1:rep, ftcs(L,T,alpha2,Nx,dt,f); end; tFTCS = toc/rep;
tic; for k=1:rep, crank_nicolson_gs(L,T,alpha2,Nx,dt,f,tol,maxIter); end; tCN = toc/rep;

Metodo_r        = {'FTCS'; 'Crank-Nicolson'};
Nx_r            = [Nx;   Nx  ];
dx_r            = [dx;   dx  ];
dt_r            = [dt;   dt  ];
r_r             = [rf;   rc  ];
err_max_r       = [errF.max_abs;       errC.max_abs      ];
err_rel_max_r   = [errF.rel_max;       errC.rel_max      ];
err_L2g_r       = [errF.L2_global;     errC.L2_global    ];
err_rL2g_r      = [errF.rel_L2_global; errC.rel_L2_global];
err_L2f_r       = [errF.L2_final;      errC.L2_final     ];
err_rL2f_r      = [errF.rel_L2_final;  errC.rel_L2_final ];
RMSE_g_r        = [errF.RMSE_global;   errC.RMSE_global  ];
RMSE_f_r        = [errF.RMSE_final;    errC.RMSE_final   ];
t_exec_r        = [tFTCS; tCN];
GS_iter_r       = [NaN; infoC.iter_prom];
obs_F = ['Expl' i_ac 'cito; estable si r <= 0.5'];
obs_C = ['Impl' i_ac 'cito; sistema resuelto con Gauss-Seidel'];
Observacion_r   = {obs_F; obs_C};

TablaResumen = table(Metodo_r, Nx_r, dx_r, dt_r, r_r, ...
    err_max_r, err_rel_max_r, ...
    err_L2g_r, err_rL2g_r, ...
    err_L2f_r, err_rL2f_r, ...
    RMSE_g_r, RMSE_f_r, ...
    t_exec_r, GS_iter_r, Observacion_r, ...
    'VariableNames', {'Metodo','Nx','dx','dt','r', ...
                      'err_max','err_rel_max', ...
                      'err_L2_global','err_rel_L2_global', ...
                      'err_L2_final','err_rel_L2_final', ...
                      'RMSE_global','RMSE_final', ...
                      't_exec_s','GS_iter_prom','Observacion'});

disp('--- Tabla resumen paper (Exp2) ---');
disp(TablaResumen);
writetable(TablaResumen, fullfile(tab_dir, 'tabla_resumen_paper.csv'));
% tabla1 como copia de la resumen para compatibilidad con entregas anteriores
writetable(TablaResumen, fullfile(tab_dir, 'tabla1_comparacion_base.csv'));

%% ===================================================================
%  EXPERIMENTO 3: ESTABILIDAD DE FTCS  (OE3)
%  Caso estable (r <= 0.5) vs caso inestable (r > 0.5).
%  Se usa malla fina (Nx=100) para que la inestabilidad sembrada por el
%  redondeo se manifieste dentro del tiempo simulado.
%  ===================================================================
Nx_e = 100;  dx_e = L / Nx_e;
r_estable   = 0.40;   dt_est = r_estable   * dx_e^2 / alpha2;
r_inestable = 0.75;   dt_ine = r_inestable * dx_e^2 / alpha2;

[xs, ts, Us, rs, ~] = ftcs(L, T, alpha2, Nx_e, dt_est, f);
[xu, tu, Uu, ru, ~] = ftcs(L, T, alpha2, Nx_e, dt_ine, f);
Ue_s = solucion_exacta(xs, ts, L, alpha2);

fig3 = figure('Color', 'w', 'Position', [100 100 900 420]);

% Subplot izquierdo: perfil final del caso estable vs exacta
subplot(1,2,1); hold on; grid on;
plot(xs, Ue_s(end,:), 'k-', 'LineWidth', 2, 'DisplayName', 'Exacta');
plot(xs, Us(end,:), 'bo--', 'MarkerSize', 4, ...
     'DisplayName', sprintf('FTCS (r = %.2f)', rs));
xlabel('x'); ylabel(sprintf('u(x, %.3f)', T));
title(sprintf('FTCS estable (r = %.2f): perfil final', rs), 'FontWeight', 'normal');
legend('Location', 'northeast', 'FontSize', 9);

% Subplot derecho: crecimiento de max|u| en escala log para el caso inestable.
% Evidencia directa de la amplificacion exponencial (modos de alta frecuencia).
subplot(1,2,2);
max_abs_Uu = max(abs(Uu), [], 2);
semilogy(tu, max_abs_Uu, 'r-', 'LineWidth', 1.4);
grid on; box on;
xlabel('Tiempo t');
ylabel('max_x |u(x,t)|  (escala log)');
title(sprintf('FTCS inestable (r = %.2f): crecimiento de max|u|', ru), ...
      'FontWeight', 'normal');

saveas(fig3, fullfile(fig_dir, 'fig3_estabilidad_ftcs.png'));

fprintf('\n[Exp3] estable  r=%.2f  max|u(T)|=%.3e (exacto %.3e)\n', ...
        rs, max(abs(Us(end,:))), max(abs(Ue_s(end,:))));
fprintf('[Exp3] inestable r=%.2f  max|u(T)|=%.3e\n', ru, max(abs(Uu(end,:))));

%% ===================================================================
%  EXPERIMENTO 4: SENSIBILIDAD / ROBUSTEZ  (OE3)
%  ===================================================================

% --- 4a) Variacion de Nx (refinamiento espacial), r fijo = 0.25 ---
% Se ajusta dt para mantener r = 0.25 constante al variar la malla espacial.
Nx_list = [25 50 100 200];
n4a     = numel(Nx_list);
rep_nx  = 3;   % repeticiones para estabilizar medicion de tiempo

Nx4_v  = zeros(n4a,1); dx4_v = zeros(n4a,1);
dt4_v  = zeros(n4a,1); r4_v  = zeros(n4a,1);
errmax_F4  = zeros(n4a,1); errL2g_F4 = zeros(n4a,1); errrel_F4 = zeros(n4a,1);
errmax_C4  = zeros(n4a,1); errL2g_C4 = zeros(n4a,1); errrel_C4 = zeros(n4a,1);
tF4 = zeros(n4a,1); tC4 = zeros(n4a,1); iterC4 = zeros(n4a,1);

for i = 1:n4a
    NxV = Nx_list(i);  dxV = L / NxV;  dtV = 0.25 * dxV^2 / alpha2;

    [xa, ta, Ua, raV, ~] = ftcs(L, T, alpha2, NxV, dtV, f);
    [~,  ~, Ub, ~, ~, ~, infoCa] = ...
        crank_nicolson_gs(L, T, alpha2, NxV, dtV, f, tol, maxIter);
    Uea = solucion_exacta(xa, ta, L, alpha2);
    eA  = calcular_error(Ua, Uea, xa, ta);
    eB  = calcular_error(Ub, Uea, xa, ta);

    tic; for kk=1:rep_nx, ftcs(L,T,alpha2,NxV,dtV,f); end
    tF4(i) = toc / rep_nx;
    tic; for kk=1:rep_nx, crank_nicolson_gs(L,T,alpha2,NxV,dtV,f,tol,maxIter); end
    tC4(i) = toc / rep_nx;

    Nx4_v(i) = NxV;  dx4_v(i) = dxV;  dt4_v(i) = dtV;  r4_v(i) = raV;
    errmax_F4(i) = eA.max_abs;    errL2g_F4(i) = eA.L2_global;  errrel_F4(i) = eA.rel_L2_global;
    errmax_C4(i) = eB.max_abs;    errL2g_C4(i) = eB.L2_global;  errrel_C4(i) = eB.rel_L2_global;
    iterC4(i) = infoCa.iter_prom;
end

% Orden de convergencia estimado: p_i = log(e_{i-1}/e_i) / log(dx_{i-1}/dx_i)
% Primera fila sin referencia anterior -> NaN.
orden_F4 = nan(n4a,1);
orden_C4 = nan(n4a,1);
for i = 2:n4a
    ratio_dx = dx4_v(i-1) / dx4_v(i);
    orden_F4(i) = log(errL2g_F4(i-1) / errL2g_F4(i)) / log(ratio_dx);
    orden_C4(i) = log(errL2g_C4(i-1) / errL2g_C4(i)) / log(ratio_dx);
end

ord_prom_F = mean(orden_F4(2:end));
ord_prom_C = mean(orden_C4(2:end));
fprintf('\n[Exp4a] Orden promedio de convergencia espacial:\n');
fprintf('        FTCS = %.3f  |  Crank-Nicolson = %.3f\n', ord_prom_F, ord_prom_C);

Tabla4a = table(Nx4_v, dx4_v, dt4_v, r4_v, ...
    errmax_F4, errL2g_F4, errrel_F4, ...
    errmax_C4, errL2g_C4, errrel_C4, ...
    orden_F4, orden_C4, tF4, tC4, iterC4, ...
    'VariableNames', {'Nx','dx','dt','r', ...
                      'err_max_FTCS','err_L2_global_FTCS','err_rel_L2_global_FTCS', ...
                      'err_max_CN',  'err_L2_global_CN',  'err_rel_L2_global_CN', ...
                      'orden_FTCS','orden_CN', ...
                      't_exec_FTCS','t_exec_CN','GS_iter_prom_CN'});

disp('--- Tabla 2 (Exp4a): sensibilidad al refinamiento espacial ---');
disp(Tabla4a);
writetable(Tabla4a, fullfile(tab_dir, 'tabla2_sensibilidad_Nx.csv'));

% Figura 4: error L2 global vs dx (log-log), convencion malla gruesa a la izquierda
dx_ref = logspace(log10(min(dx4_v)*0.7), log10(max(dx4_v)*1.3), 30);
ref_p2 = errL2g_F4(1) * (dx_ref / dx4_v(1)).^2;

fig4 = figure('Color', 'w', 'Position', [100 100 720 480]);
loglog(dx4_v, errL2g_F4, 'bo-', 'LineWidth', 1.4, 'DisplayName', 'FTCS');
hold on; grid on;
loglog(dx4_v, errL2g_C4, 'rs-', 'LineWidth', 1.4, 'DisplayName', 'Crank-Nicolson');
loglog(dx_ref, ref_p2, 'k--', 'LineWidth', 0.8, 'DisplayName', 'Referencia p = 2');
set(gca, 'XDir', 'reverse');   % malla gruesa (dx grande) a la izquierda
xlabel('Paso espacial dx  (dx = L/Nx)'); ylabel('Error L2 global');
title('Convergencia espacial: error L2 global vs dx', 'FontWeight', 'normal');
legend('Location', 'northwest', 'FontSize', 9);
saveas(fig4, fullfile(fig_dir, 'fig4_convergencia_Nx.png'));

% --- 4b) Variacion de r (dt), Nx fijo = 50 ---
% Incluye r = 0.75 (inestable) para conectar con fig3_estabilidad_ftcs.png.
% Para FTCS inestable: errores se omiten (NaN) y se reporta max|u| en su lugar.
r_list_b = [0.10 0.25 0.40 0.50 0.75];
n4b = numel(r_list_b);

r4b_v    = zeros(n4b,1); dt4b_v  = zeros(n4b,1);
estado_F = cell(n4b,1);
emaxF4b  = nan(n4b,1);  eL2gF4b = nan(n4b,1);  erelF4b = nan(n4b,1);
maxabsF4b = nan(n4b,1);
emaxC4b  = nan(n4b,1);  eL2gC4b = nan(n4b,1);  erelC4b = nan(n4b,1);
iterC4b  = nan(n4b,1);
obs4b    = cell(n4b,1);

for i = 1:n4b
    rV  = r_list_b(i);
    dtV = rV * dx^2 / alpha2;
    [xa, ta, Ua, raV, ~] = ftcs(L, T, alpha2, Nx, dtV, f);
    [~, ~, Ub, ~, ~, ~, infoB] = ...
        crank_nicolson_gs(L, T, alpha2, Nx, dtV, f, tol, maxIter);
    Uea = solucion_exacta(xa, ta, L, alpha2);
    eB  = calcular_error(Ub, Uea, xa, ta);

    r4b_v(i)  = raV;   dt4b_v(i) = dtV;
    emaxC4b(i)  = eB.max_abs;
    eL2gC4b(i)  = eB.L2_global;
    erelC4b(i)  = eB.rel_L2_global;
    iterC4b(i)  = infoB.iter_prom;

    if rV <= 0.5
        eA = calcular_error(Ua, Uea, xa, ta);
        estado_F{i} = 'estable';
        emaxF4b(i)  = eA.max_abs;
        eL2gF4b(i)  = eA.L2_global;
        erelF4b(i)  = eA.rel_L2_global;
        obs4b{i}    = '';
    else
        % r > 0.5: FTCS viola el criterio de Von Neumann (r <= 0.5).
        % Los errores no se interpretan como metrica de precision; se registra
        % max|u| para evidenciar la explosion numerica.
        estado_F{i}  = 'inestable';
        maxabsF4b(i) = max(abs(Ua(end,:)));
        obs4b{i}     = 'Viola r<=0.5; max|u| refleja crecimiento numerico';
    end
end

Tabla4b = table(r4b_v, dt4b_v, estado_F, ...
    emaxF4b, eL2gF4b, erelF4b, maxabsF4b, ...
    emaxC4b, eL2gC4b, erelC4b, iterC4b, obs4b, ...
    'VariableNames', {'r','dt','estado_FTCS', ...
                      'err_max_FTCS','err_L2_global_FTCS','err_rel_L2_global_FTCS', ...
                      'max_abs_u_FTCS', ...
                      'err_max_CN','err_L2_global_CN','err_rel_L2_global_CN', ...
                      'CN_iter_prom','observacion'});

disp('--- Tabla 3 (Exp4b): sensibilidad al paso temporal (r) ---');
disp(Tabla4b);
writetable(Tabla4b, fullfile(tab_dir, 'tabla3_sensibilidad_r.csv'));

% --- 4c) Variacion de alpha2, Nx=50, r fijo = 0.25 ---
% Se ajusta dt para mantener r constante y analizar la difusividad bajo
% una condicion de estabilidad uniforme: dt = 0.25 * dx^2 / alpha2.
a2_list = [0.005 0.01 0.05 0.10];
n4c     = numel(a2_list);
rep_a2  = 3;

a2_v  = zeros(n4c,1); dt4c_v = zeros(n4c,1); r4c_v  = zeros(n4c,1);
emaxF4c = zeros(n4c,1); eL2gF4c = zeros(n4c,1); erelF4c = zeros(n4c,1);
emaxC4c = zeros(n4c,1); eL2gC4c = zeros(n4c,1); erelC4c = zeros(n4c,1);
tF4c = zeros(n4c,1); tC4c = zeros(n4c,1); iterC4c = zeros(n4c,1);

for i = 1:n4c
    a2V = a2_list(i);  dtV = 0.25 * dx^2 / a2V;
    [xa, ta, Ua, raV, ~] = ftcs(L, T, a2V, Nx, dtV, f);
    [~, ~, Ub, ~, ~, ~, infoCc] = ...
        crank_nicolson_gs(L, T, a2V, Nx, dtV, f, tol, maxIter);
    Uea = solucion_exacta(xa, ta, L, a2V);
    eA  = calcular_error(Ua, Uea, xa, ta);
    eB  = calcular_error(Ub, Uea, xa, ta);

    tic; for kk=1:rep_a2, ftcs(L,T,a2V,Nx,dtV,f); end
    tF4c(i) = toc / rep_a2;
    tic; for kk=1:rep_a2, crank_nicolson_gs(L,T,a2V,Nx,dtV,f,tol,maxIter); end
    tC4c(i) = toc / rep_a2;

    a2_v(i) = a2V;   dt4c_v(i) = dtV;   r4c_v(i) = raV;
    emaxF4c(i) = eA.max_abs;  eL2gF4c(i) = eA.L2_global;  erelF4c(i) = eA.rel_L2_global;
    emaxC4c(i) = eB.max_abs;  eL2gC4c(i) = eB.L2_global;  erelC4c(i) = eB.rel_L2_global;
    iterC4c(i) = infoCc.iter_prom;
end

Tabla4c = table(a2_v, dt4c_v, r4c_v, ...
    emaxF4c, eL2gF4c, erelF4c, ...
    emaxC4c, eL2gC4c, erelC4c, ...
    tF4c, tC4c, iterC4c, ...
    'VariableNames', {'alpha2','dt','r', ...
                      'err_max_FTCS','err_L2_global_FTCS','err_rel_L2_global_FTCS', ...
                      'err_max_CN',  'err_L2_global_CN',  'err_rel_L2_global_CN', ...
                      't_exec_FTCS','t_exec_CN','GS_iter_prom_CN'});

disp('--- Tabla 4 (Exp4c): sensibilidad a la difusividad alpha2 ---');
disp(Tabla4c);
writetable(Tabla4c, fullfile(tab_dir, 'tabla4_sensibilidad_alpha2.csv'));

fprintf('\nListo. Figuras en %s\nTablas en %s\n', fig_dir, tab_dir);
