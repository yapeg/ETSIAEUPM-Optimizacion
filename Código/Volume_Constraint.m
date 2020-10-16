function [VC_ineq, VC_eq] = Volume_Constraint(F, z, F0, F1, V_FA)
    % Se obliga a que el volumen sea el mismo que el proporcionado por la
    % solución analítica (restricción estricta).

    F = [F0 F F1];

    VC_ineq = [];                        % Restricción unilateral (no existe)

    I_V = F.^2;

    % Diferencia entre el volumen introducido y el obtenido por la regla del
    % trapecio.
    VC_eq = trapz(z, I_V) - V_FA;        % Restricción estricta
    
end
