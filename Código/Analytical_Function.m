%% Soluci�n anal�tica del problema

function [FA, A_FA, V_FA] = Analytical_Function(z, F0, F1)
    % F0 y F1 sirven como condiciones de contorno

    N = max(size(z)) - 1;

    % Se resuelve el sistema de inc�gnitas con las condiciones de contorno.
    CC = fsolve(@(CC) [CC(1) * cosh((z(1) - CC(2)) / CC(1)) - F0, CC(1) * cosh((z(N+1) - CC(2)) / CC(1)) - F1], [1, 1]);
    
    % Array funci�n anal�tica.
    FA = CC(1) * cosh((z - CC(2)) / CC(1));
    
    % C�lculo del �rea y volumen con la funci�n anal�tica
    syms x

    F = CC(1) * cosh((x - CC(2)) / CC(1));

    % Integral del �rea
    I_A = F * sqrt(1 + (diff(F))^2);
    A_FA = double(int(I_A, [z(1) z(N+1)]));

    % Integral del volumen
    I_V = F^2;
    V_FA = double(int(I_V, [z(1) z(N+1)]));

end



