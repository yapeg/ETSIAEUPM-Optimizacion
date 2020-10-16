function EULER = Euler_Equation(F, z, F0, F1)

    N = max(size(z)) - 1;
    F = [F0 F F1];

    d2F = zeros(1, N+1);
    
    % Diferencias finitas centradas, para derivadas de primer y segundo
    % orden.
    
    for i = 1:N+1
        % extremo inferior obliga a utilizar DFs adelantadas.        
        if i == 1
            dF(i) = (F(i+1) - F(i)) / (z(i+1) - z(i));
            d2F(i) = (F(i+2) - 2 * F(i+1) + F(i)) / ((z(i+2) - z(i+1)) * (z(i+1) - z(i)));
        % extremo superior obliga a utilizar DFs atrasadas.
        elseif i == N+1
            dF(i) = (F(i) - F(i-1)) / (z(i) - z(i-1));
            d2F(i) = (F(i) - 2 * F(i-1) + F(i-2)) / ((z(i) - z(i-1)) * (z(i-1) - z(i-2)));
        
        else
            dF(i) = (F(i+1) - F(i-1)) / (z(i+1) - z(i-1));
            d2F(i) = (F(i+1) - 2 * F(i) + F(i-1)) / ((z(i+1) - z(i)) * (z(i) - z(i-1)));
        end
      
    end
    
    EULER = F .* d2F - (1 + dF.^2);
end