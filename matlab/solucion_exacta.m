function U_exacta = solucion_exacta(x, t, L, alpha2)
% SOLUCION_EXACTA  Solucion analitica de la ecuacion del calor 1D.
%
%   Proposito:
%       Evaluar la solucion exacta de la ecuacion del calor
%           du/dt = alpha2 * d2u/dx2,  u(0,t)=u(L,t)=0,  u(x,0)=sin(pi*x/L)
%       cuya forma cerrada es
%           u(x,t) = exp(-alpha2*(pi/L)^2 * t) * sin(pi*x/L).
%       Se usa como referencia para medir el error numerico de FTCS y
%       Crank-Nicolson nodo por nodo y tiempo por tiempo.
%
%   Entradas:
%       x       - vector fila de posiciones espaciales (1 x (Nx+1)).
%       t       - vector de instantes temporales (1 x (Nt+1)).
%       L       - longitud de la barra.
%       alpha2  - difusividad termica (alpha^2).
%
%   Salidas:
%       U_exacta - matriz (Nt+1) x (Nx+1). Cada fila es el perfil espacial
%                  en un instante t(n); cada columna la evolucion temporal
%                  en una posicion x(i). Misma orientacion que las matrices
%                  U devueltas por ftcs.m y crank_nicolson_gs.m.
%
%   Justificacion numerica:
%       Disponer de la solucion analitica permite cuantificar de forma
%       objetiva el error de discretizacion y evaluar convergencia, en
%       lugar de comparar esquemas solo entre si.

    % --- Validacion de entradas ---
    if L <= 0
        error('solucion_exacta:L', 'L debe ser positivo.');
    end
    if alpha2 < 0
        error('solucion_exacta:alpha2', 'alpha2 no puede ser negativo.');
    end

    % Se fuerza orientacion: t como columna, x como fila -> producto externo.
    t = t(:);          % (Nt+1) x 1
    x = x(:).';        % 1 x (Nx+1)

    decaimiento = exp(-alpha2 * (pi / L)^2 * t);   % (Nt+1) x 1
    perfil      = sin(pi * x / L);                 % 1 x (Nx+1)

    U_exacta = decaimiento * perfil;               % (Nt+1) x (Nx+1)
end
