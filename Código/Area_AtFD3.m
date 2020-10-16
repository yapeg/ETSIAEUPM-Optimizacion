function A = Area_AtFD3(F, z, F0, F1)
    % Se busca optimizar el área.
    
    N = max(size(z)) - 1;
    
    F = [F0 F F1];
    dF = zeros(1, N+1);
    
    % Diferencias finitas atrasadas, con tres puntos.
    
    for i = 1:N+1
        % extremo inferior obliga a utilizar derivada adelantada.
        if i == 1
            dF(i) = (F(i+1) - F(i)) / (z(i+1) - z(i));
        % segundo extremo inferior obliga a utilizar derivada centrada.
        elseif i == 2
            dF(i) = (F(i+1) - F(i-1)) / (z(i+1) - z(i-1));
            
        else
            dF(i) = (3 * F(i) - 4 * F(i-1) + F(i-2)) / (z(i) - z(i-2));
        end
        
    end
    
    % Integrando
    I = F .* sqrt(1 + dF.^2);
    
    % Integral por método del trapecio.
    A = trapz(z, I);
    
end
    
    