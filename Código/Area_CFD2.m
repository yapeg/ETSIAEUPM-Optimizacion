function A = Area_CFD2(F, z, F0, F1)
    % Se busca optimizar el área.

    N = max(size(z)) - 1;

    F = [F0 F F1];
    dF = zeros(1, N+1);

    % Diferencias finitas centradas, con dos puntos.
    for i = 1:N+1
        % extremo inferior obliga a utilizar derivada adelantada.
        if i == 1
            dF(i) = (F(i+1) - F(i)) / (z(i+1) - z(i));
        % extremo superior obliga a utilizar derivada atrasada.
        elseif i == N+1
            dF(i) = (F(i) - F(i-1)) / (z(i) - z(i-1));
        
        else    
            dF(i) = (F(i+1) - F(i-1)) / (z(i+1) - z(i-1));
        end
    
    end

    % Integrando.
    I = F .* sqrt(1 + dF.^2);

    % Integral por método del trapecio.
    A = trapz(z, I);

end