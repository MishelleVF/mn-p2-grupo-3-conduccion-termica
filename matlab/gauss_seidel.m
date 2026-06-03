function [x, iter, err, converge] = gauss_seidel(A, b, x0, tol, maxIter)
% GAUSS_SEIDEL  Metodo iterativo de Gauss-Seidel para A*x = b.
%
%   Proposito:
%       Resolver de forma iterativa el sistema lineal A*x = b. En este
%       proyecto se usa para resolver el sistema que aparece en CADA paso
%       temporal del esquema implicito de Crank-Nicolson.
%
%   Entradas:
%       A       - matriz de coeficientes (n x n), con diagonal no nula.
%       b       - vector de terminos independientes (n x 1).
%       x0      - aproximacion inicial (n x 1).
%       tol     - tolerancia para el criterio de parada.
%       maxIter - numero maximo de iteraciones.
%
%   Salidas:
%       x        - solucion aproximada (n x 1).
%       iter     - numero de iteraciones realizadas.
%       err      - error relativo final (norma infinito).
%       converge - true si se alcanzo la tolerancia antes de maxIter.
%
%   Justificacion numerica:
%       Gauss-Seidel actualiza cada componente usando los valores ya
%       calculados en la misma iteracion:
%           x_i^{(k+1)} = (1/a_ii) [ b_i - sum_{j<i} a_ij x_j^{(k+1)}
%                                         - sum_{j>i} a_ij x_j^{(k)} ].
%       La matriz de Crank-Nicolson es tridiagonal con diagonal (1+r) y
%       fuera de la diagonal (-r/2), por lo que es estrictamente diagonal
%       dominante (1+r > r) y simetrica definida positiva: Gauss-Seidel
%       converge. El criterio de parada usa la norma infinito relativa
%       del cambio entre iteraciones, robusta frente a la escala.

    % --- Validacion de dimensiones ---
    [m, n] = size(A);
    if m ~= n
        error('gauss_seidel:cuadrada', 'A debe ser cuadrada.');
    end
    b  = b(:);
    x0 = x0(:);
    if numel(b) ~= n || numel(x0) ~= n
        error('gauss_seidel:dims', 'b y x0 deben tener longitud %d.', n);
    end
    if any(abs(diag(A)) < eps)
        error('gauss_seidel:diag', 'A tiene un elemento diagonal nulo.');
    end

    x        = x0;
    converge = false;
    err      = Inf;

    for iter = 1:maxIter
        x_old = x;
        for i = 1:n
            s1 = A(i, 1:i-1)   * x(1:i-1);       % componentes ya actualizadas
            s2 = A(i, i+1:n)   * x_old(i+1:n);   % componentes de la iter previa
            x(i) = (b(i) - s1 - s2) / A(i, i);
        end
        % Error relativo en norma infinito
        err = norm(x - x_old, inf) / max(norm(x, inf), eps);
        if err < tol
            converge = true;
            break;
        end
    end
end
