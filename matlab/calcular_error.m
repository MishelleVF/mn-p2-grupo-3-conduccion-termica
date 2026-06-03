function err = calcular_error(U_num, U_exacta, x, t)
% CALCULAR_ERROR  Metricas de error de una solucion numerica frente a la exacta.
%
%   Proposito:
%       Cuantificar la precision de FTCS o Crank-Nicolson comparando su
%       matriz U_num con la solucion exacta U_exacta, ambas con la misma
%       orientacion (Nt+1) x (Nx+1).
%
%   Entradas:
%       U_num    - solucion numerica (Nt+1) x (Nx+1).
%       U_exacta - solucion exacta   (Nt+1) x (Nx+1).
%       x        - nodos espaciales (para el peso dx en las normas L2).
%       t        - instantes temporales (para el peso dt y el error por tiempo).
%
%   Salidas (struct err):
%       err.max_abs        - max|E|, error puntual maximo en todo el dominio.
%       err.L2_global      - norma L2 espacio-tiempo: sqrt(dx*dt*sum(E^2)).
%       err.L2_final       - norma L2 en t=T:         sqrt(dx*sum(E(end,:)^2)).
%       err.rel_max        - max_abs / max|U_exacta|.
%       err.rel_L2_global  - L2_global / ||U_exacta||_L2_global.
%       err.rel_L2_final   - L2_final  / ||U_exacta(T)||_L2.
%       err.RMSE_global    - RMSE sobre todos los nodos espacio-temporales.
%       err.RMSE_final     - RMSE sobre los nodos espaciales en t=T.
%       err.L2_por_tiempo  - (Nt+1)x1: norma L2 espacial en cada instante.
%
%   Justificacion numerica:
%       Se combinan metricas L-infinito, L2 y RMSE ponderadas por dx y dt
%       para que sean comparables entre mallas distintas. Los errores
%       relativos normalizan respecto a la magnitud de la solucion exacta,
%       facilitando la comparacion entre casos con distintos alpha2.

    % --- Validacion de dimensiones ---
    if ~isequal(size(U_num), size(U_exacta))
        error('calcular_error:dims', ...
              'U_num y U_exacta deben tener el mismo tamano.');
    end

    E  = U_num - U_exacta;            % matriz de errores
    dx = x(2) - x(1);
    if numel(t) > 1
        dt = t(2) - t(1);
    else
        dt = 1;
    end

    % --- Metricas absolutas ---
    err.max_abs   = max(abs(E(:)));
    err.L2_global = sqrt(dx * dt * sum(E(:).^2));
    Ef = E(end, :);
    err.L2_final  = sqrt(dx * sum(Ef.^2));

    % --- Error relativo maximo ---
    denom = max(abs(U_exacta(:)));
    if denom < eps; denom = 1; end
    err.rel_max = err.max_abs / denom;

    % --- Error relativo L2 global (normalizado por la norma de la exacta) ---
    denom_g = sqrt(dx * dt * sum(U_exacta(:).^2));
    if denom_g < eps; denom_g = 1; end
    err.rel_L2_global = err.L2_global / denom_g;

    % --- Error relativo L2 en el tiempo final ---
    denom_f = sqrt(dx * sum(U_exacta(end,:).^2));
    if denom_f < eps; denom_f = 1; end
    err.rel_L2_final = err.L2_final / denom_f;

    % --- RMSE global (promedio cuadratico sobre todos los nodos espacio-tiempo) ---
    err.RMSE_global = sqrt(sum(E(:).^2) / numel(E));

    % --- RMSE en el tiempo final ---
    err.RMSE_final = sqrt(sum(Ef.^2) / numel(Ef));

    % --- Norma L2 espacial por cada instante temporal ---
    err.L2_por_tiempo = sqrt(dx * sum(E.^2, 2));   % (Nt+1) x 1
end
