start = tic;
%% Autores:
% Fernández Eguizábal, Andrés
% Galván Castro, Adal David
% Palmeiro Sánchez, Sergio
% Pego Martínez, Yago

%% FORMA DE EQUILIBRIO DE LA SUPERFICIE LIBRE DE UN FLUIDO
% Optimización

% Programa principal

%% GENERALIDADES
    N = 50;                     % nº de puntos
    delta = 1/N;

    z_0 = 0;                    % z inicial
    z_f = 1;                    % z final
    z = z_0:delta:z_f;          % nube de puntos (N+1)

    count = 1;                  % contador

% Condiciones del problema
    r0 = 1;                     % radio de los anillos
    F = ones(1, N-1);           % radios intermedios (inicialización)
    F0 = r0;                    % radio inferior
    F1 = r0;                    % radio superior

%% 1. FUNCIÓN ANALÍTICA

% Se obtiene la función FA, el área A_FA y volumen v_FA de forma analítica
% del problema presentado, con las condiciones dadas.
    [FA, A_FA, V_FA] = Analytical_Function(z, F0, F1);

% Ploteado de la solución analítica. Variación del radio con la altura.
    figure(count)
    plot(FA, z)
    axis([0.7 1, 0 1])
    title('Solución analítica')
    xlabel('r')
    ylabel('z')
    grid on
    
    
%% 2. OPTIMIZACIÓN CON MÉTODO DFP

% Método de optimización basado en el gradiente.
% Desarrollado por Davidon, Fletcher y Powell, en 1963.

% Opciones de optimización para la búsqueda de mínimos mediante DFP.
    options = optimoptions('fminunc');
    options = optimoptions(options, 'Display', 'off');
    %options = optimoptions(options, 'Display', 'iter');
    options = optimoptions(options, 'MaxFunctionEvaluations', 3e5);           
    options = optimoptions(options, 'MaxIterations', 5e4);               
    options = optimoptions(options, 'FunctionTolerance', 1e-8);                 
    options = optimoptions(options, 'OptimalityTolerance', 1e-8);
    options = optimoptions(options, 'StepTolerance', 1e-8);
    options = optimoptions(options, 'Algorithm', 'quasi-newton');    
    options = optimoptions(options, 'HessUpdate', 'DFP');
    
% Optimización DFP
    tic;
    [F_DFP_CFD, A_DFP_CFD] = fminunc(@(F) Area_CFD2(F, z, F0, F1), F, options);
    t_DFP_CFD = toc;
    
    % Cálculo del volumen DFP
    F_DFP_CFD = [F0 F_DFP_CFD F1];
    V_DFP_CFD = trapz(z, F_DFP_CFD.^2);
   
% Ploteado del resultado obtenido mediante optimización DFP.
    count = count + 1;
    figure(count)
    plot(F_DFP_CFD, z);
    axis([0.7 1, 0 1])
    title('Optimización por método DFP')
    xlabel('r')
    ylabel('z')
    grid on
    
% Cálculo del error del método DFP, DFs centradas.
    E_DFP_CFD = abs(F_DFP_CFD - FA);
    count = count + 1;
    figure(count)
    plot(E_DFP_CFD, z);
    title('Error del método DFP')
    xlabel('Error DFP')
    ylabel('z')
    grid on
    
 % Variación de los parámetros del problema.
 
    % Cambio en los radios de los anillos (condiciones de contorno).
    % Vector de valores.
    r0_anillos = [0.8 0.9 1 1.25 1.5];
        
        % Para F0 = 0.8:
        F01 = ['F_0 = ', num2str(r0_anillos(1))];
        F0_CC = r0_anillos(1);
        F1_CC = r0_anillos(1);
        F_CC(1:N-1) = r0_anillos(1);
        tic;
        [F_DFP_temp, A_DFP_CC(1)] = fminunc(@(F) Area_CFD2(F, z, F0_CC, F1_CC), F_CC, options);
        t_DFP_CC(1) = toc;
        F_DFP_CC(1,:) = [F0_CC F_DFP_temp F1_CC];
           
        % Para F0 = 0.9:
        F02 = ['F_0 = ', num2str(r0_anillos(2))];
        F0_CC = r0_anillos(2);
        F1_CC = r0_anillos(2);
        F_CC(1:N-1) = r0_anillos(2);
        tic;
        [F_DFP_temp, A_DFP_CC(2)] = fminunc(@(F) Area_CFD2(F, z, F0_CC, F1_CC), F_CC, options); 
        t_DFP_CC(2) = toc;
        F_DFP_CC(2,:) = [F0_CC F_DFP_temp F1_CC];   

        % Para F0 = 1:
        F03 = ['F_0 = ', num2str(r0_anillos(3))];
        F0_CC = r0_anillos(3);
        F1_CC = r0_anillos(3);
        F_CC(1:N-1) = r0_anillos(3);
        tic;
        [F_DFP_temp, A_DFP_CC(3)] = fminunc(@(F) Area_CFD2(F, z, F0_CC, F1_CC), F_CC, options); 
        t_DFP_CC(3) = toc;
        F_DFP_CC(3,:) = [F0_CC F_DFP_temp F1_CC];
        
        % Para F0 = 1.25:
        F04 = ['F_0 = ', num2str(r0_anillos(4))];
        F0_CC = r0_anillos(4);
        F1_CC = r0_anillos(4);
        F_CC(1:N-1) = r0_anillos(4);
        tic;
        [F_DFP_temp, A_DFP_CC(4)] = fminunc(@(F) Area_CFD2(F, z, F0_CC, F1_CC), F_CC, options); 
        t_DFP_CC(4) = toc;
        F_DFP_CC(4,:) = [F0_CC F_DFP_temp F1_CC];
        
        % Para F0 = 1.5:
        F05 = ['F_0 = ', num2str(r0_anillos(5))];
        F0_CC = r0_anillos(5);
        F1_CC = r0_anillos(5);
        F_CC(1:N-1) = r0_anillos(5);
        tic;
        [F_DFP_temp, A_DFP_CC(5)] = fminunc(@(F) Area_CFD2(F, z, F0_CC, F1_CC), F_CC, options); 
        t_DFP_CC(5) = toc;
        F_DFP_CC(5,:) = [F0_CC F_DFP_temp F1_CC];
        
    % Comparación para las distintas condiciones de contorno.
    count = count + 1;
    figure(count)
    plot(F_DFP_CC(1,:), z, F_DFP_CC(2,:), z, F_DFP_CC(3,:), z, F_DFP_CC(4,:), z, F_DFP_CC(5,:), z);
    axis([0 1.5 0 1])
    title('Influencia de las condiciones de contorno (DFP)')
    legend(F01, F02, F03, F04, F05)
    xlabel('r')
    ylabel('z')
    grid on

    
    % Cambio en el número de puntos.
    M = [10, 50, 100, 250];
    
        % Para 10 puntos:
        M1 = ['N = ', num2str(M(1))];
        dz = 1/M(1);
        z_10 = 0:dz:1;
        FM = ones(1, M(1) - 1);
        tic;
        [F_DFP_10, A_DFP_M(1)] = fminunc(@(F) Area_CFD2(F, z_10, F0, F1), FM, options);
        t_DFP_M(1) = toc;
        F_DFP_10 = [F0 F_DFP_10 F1];
        
        % Para 50 puntos:
        M2 = ['N = ', num2str(M(2))];
        dz = 1/M(2);
        z_50 = 0:dz:1;
        FM = ones(1, M(2) - 1);
        tic;
        [F_DFP_50, A_DFP_M(2)] = fminunc(@(F) Area_CFD2(F, z_50, F0, F1), FM, options);
        t_DFP_M(2) = toc;
        F_DFP_50 = [F0 F_DFP_50 F1];        

        % Para 100 puntos:
        M3 = ['N = ', num2str(M(3))];
        dz = 1/M(3);
        z_100 = 0:dz:1;
        FM = ones(1, M(3) - 1);
        tic;
        [F_DFP_100, A_DFP_M(3)] = fminunc(@(F) Area_CFD2(F, z_100, F0, F1), FM, options);
        t_DFP_M(3) = toc;
        F_DFP_100 = [F0 F_DFP_100 F1];

        % Para 250 puntos:
        M4 = ['N = ', num2str(M(4))];
        dz = 1/M(4);
        z_250 = 0:dz:1;
        FM = ones(1, M(4) - 1);
        tic;
        [F_DFP_250, A_DFP_M(4)] = fminunc(@(F) Area_CFD2(F, z_250, F0, F1), FM, options);
        t_DFP_M(4) = toc;
        F_DFP_250 = [F0 F_DFP_250 F1];
        
            
    % Comparación para distinto número de puntos.
    count = count + 1;
    figure(count)
    plot(F_DFP_10, z_10, F_DFP_50, z_50, F_DFP_100, z_100, F_DFP_250, z_250);
    axis([0.8 1, 0 1])
    title('Influencia del número de puntos (DFP)')
    legend(M1, M2, M3, M4)
    xlabel('r')
    ylabel('z')
    grid on

        
    % Cambio en la distribución de los puntos.
    % Distribución sinusoidal simple.
    z_sin = sin((0:delta:1) * pi / 2);
        
    % Optimización DFP para distribución sinusoidal de puntos.
    tic;
    [F_DFP_sin, A_DFP_sin] = fminunc(@(F) Area_CFD2(F, z_sin, F0, F1), F, options);
    t_DFP_sin = toc;
    F_DFP_sin = [F0 F_DFP_sin F1];
       
        % Ploteado para la distribución sinusoidal.
        count = count + 1;
        figure(count)
        plot(F_DFP_sin, z_sin);
        axis([0.7 1, 0 1])
        title('Optimización por método DFP: distribución sinusoidal de puntos')
        xlabel('r')
        ylabel('z')
        grid on
       
        % Error de la distribución sinusoidal.
        % Obtención de la solución analítica para la distribución sinusoidal.
        [FA_sin, A_FA_sin, V_FA_sin] = Analytical_Function(z_sin, F0, F1);
        
        E_DFP_sin = abs(F_DFP_sin - FA_sin);
        count = count + 1;
        figure(count)
        plot(E_DFP_sin, z_sin);
        title('Error del método DFP: distribución sinusoidal')
        xlabel('Error DFP')
        ylabel('z')
        grid on
    
    % Comparación de las distintas distribuciones de puntos.
    count = count + 1;
    figure(count)
    plot(F_DFP_CFD, z, F_DFP_sin, z_sin);
    axis([0.8 1, 0 1])
    title('Influencia de la distribución de los puntos (DFP)')
    legend('Distribución equiespaciada', 'Distribución sinusoidal')
    xlabel('r')
    ylabel('z')
    grid on
        

    % Cambio en la discretización de las derivadas.
    tic;
    [F_DFP_AtFD, A_DFP_AtFD] = fminunc(@(F) Area_AtFD3(F, z, F0, F1), F, options);
    t_DFP_AtFD = toc;

        % Ploteado del resultado obtenido mediante optimización DFP.
        count = count + 1;
        figure(count)
        F_DFP_AtFD = [F0 F_DFP_AtFD F1];
        plot(F_DFP_AtFD, z);
        axis([0.7 1, 0 1])
        title('Optimización por método DFP: DFs atrasadas')
        xlabel('r')
        ylabel('z')
        grid on
    
        % Cálculo del error del método DFP, DFs atrasadas.
        E_DFP_AtFD = abs(F_DFP_AtFD - FA);
    
        count = count + 1;
        figure(count)
        plot(E_DFP_AtFD, z);
        title('Error del método DFP: DFs atrasadas')
        xlabel('Error DFP')
        ylabel('z')
        grid on

    % Comparación entre discretizaciones.
    count = count + 1;
    figure(count)
    plot(F_DFP_CFD, z, F_DFP_AtFD, z);
    axis([0.9 1, 0 0.1])
    title('Influencia de la discretización (DFP)')
    legend('DFs Centradas', 'DFs Atrasadas')
    xlabel('r')
    ylabel('z')
    grid on    

    
    
%% 3. OPTIMIZACIÓN CON MÉTODO BFGS

% Método de optimización basado en el gradiente.
% Desarrollado por Broyden, Fletcher, Goldfarb y Shanno, en 1970.

% Opciones de optimización para la búsqueda de mínimos mediante BFGS.
    options = optimoptions('fminunc');
    options = optimoptions(options, 'Display', 'off');
    %options = optimoptions(options, 'Display', 'iter');
    options = optimoptions(options, 'MaxFunctionEvaluations', 2e6);
    options = optimoptions(options, 'MaxIterations', 1e4);
    options = optimoptions(options, 'FunctionTolerance', 1e-8);
    options = optimoptions(options, 'OptimalityTolerance', 1e-8);
    options = optimoptions(options, 'StepTolerance', 1e-8);
    options = optimoptions(options, 'Algorithm', 'quasi-newton');
    options = optimoptions(options, 'ObjectiveLimit', 1e-6);
    options = optimoptions(options, 'HessUpdate', 'BFGS');
    
% Optimización BFGS 
    tic;
    [F_BFGS_CFD, A_BFGS_CFD] = fminunc(@(F) Area_CFD2(F, z, F0, F1), F, options);
    t_BFGS_CFD = toc;
    
    % Cálculo del volumen BFGS.
    F_BFGS_CFD = [F0 F_BFGS_CFD F1];
    V_BFGS_CFD = trapz(z, F_BFGS_CFD.^2);
    
% Ploteado del resultado obtenido mediante optimización BFGS.
    count = count + 1;
    figure(count)
    plot(F_BFGS_CFD, z);
    axis([0.7 1, 0 1])
    title('Optimización por método BFGS')
    xlabel('r')
    ylabel('z')
    grid on

% Cálculo del error del método BFGS, DFs centradas.
    E_BFGS_CFD = abs(F_BFGS_CFD - FA);
    count = count + 1;
    figure(count)
    plot(E_BFGS_CFD, z);
    title('Error del método BFGS')
    xlabel('Error BFGS')
    ylabel('z')
    grid on

% Variación de los parámetros del problema.

   % Cambio en los radios de los anillos (condiciones de contorno).
   % Vector de valores.
   r0_anillos = [0.8 0.9 1 1.25 1.5];
        
        % Para F0 = 0.8:
        F01 = ['F_0 = ', num2str(r0_anillos(1))];
        F0_CC = r0_anillos(1);
        F1_CC = r0_anillos(1);
        F_CC(1:N-1) = r0_anillos(1);
        tic;
        [F_BFGS_temp, A_BFGS_CC(1)] = fminunc(@(F) Area_CFD2(F, z, F0_CC, F1_CC), F_CC, options); 
        t_BFGS_CC(1) = toc;
        F_BFGS_CC(1,:) = [F0_CC F_BFGS_temp F1_CC];
           
        % Para F0 = 0.9:
        F02 = ['F_0 = ', num2str(r0_anillos(2))];
        F0_CC = r0_anillos(2);
        F1_CC = r0_anillos(2);
        F_CC(1:N-1) = r0_anillos(2);
        tic;
        [F_BFGS_temp, A_BFGS_CC(2)] = fminunc(@(F) Area_CFD2(F, z, F0_CC, F1_CC), F_CC, options); 
        t_BFGS_CC(2) = toc;
        F_BFGS_CC(2,:) = [F0_CC F_BFGS_temp F1_CC];   

        % Para F0 = 1:
        F03 = ['F_0 = ', num2str(r0_anillos(3))];
        F0_CC = r0_anillos(3);
        F1_CC = r0_anillos(3);
        F_CC(1:N-1) = r0_anillos(3);
        tic;
        [F_BFGS_temp, A_BFGS_CC(3)] = fminunc(@(F) Area_CFD2(F, z, F0_CC, F1_CC), F_CC, options); 
        t_BFGS_CC(3) = toc;
        F_BFGS_CC(3,:) = [F0_CC F_BFGS_temp F1_CC];
        
        % Para F0 = 1.25:
        F04 = ['F_0 = ', num2str(r0_anillos(4))];
        F0_CC = r0_anillos(4);
        F1_CC = r0_anillos(4);
        F_CC(1:N-1) = r0_anillos(4);
        tic;
        [F_BFGS_temp, A_BFGS_CC(4)] = fminunc(@(F) Area_CFD2(F, z, F0_CC, F1_CC), F_CC, options); 
        t_BFGS_CC(4) = toc;
        F_BFGS_CC(4,:) = [F0_CC F_BFGS_temp F1_CC];
        
        % Para F0 = 1.5:
        F05 = ['F_0 = ', num2str(r0_anillos(5))];
        F0_CC = r0_anillos(5);
        F1_CC = r0_anillos(5);
        F_CC(1:N-1) = r0_anillos(5);
        tic;
        [F_BFGS_temp, A_BFGS_CC(5)] = fminunc(@(F) Area_CFD2(F, z, F0_CC, F1_CC), F_CC, options); 
        t_BFGS_CC(5) = toc;
        F_BFGS_CC(5,:) = [F0_CC F_BFGS_temp F1_CC];
        
    % Comparación para las distintas condiciones de contorno.
    count = count + 1;
    figure(count)
    plot(F_BFGS_CC(1,:), z, F_BFGS_CC(2,:), z, F_BFGS_CC(3,:), z, F_BFGS_CC(4,:), z, F_BFGS_CC(5,:), z);
    axis([0 1.5 0 1])
    title('Influencia de las condiciones de contorno (BFGS)')
    legend(F01, F02, F03, F04, F05)
    xlabel('r')
    ylabel('z')
    grid on
    
    
    % Cambio en el número de puntos.
    M = [10, 50, 100, 250];
    
        % Para 10 puntos:
        M1 = ['N = ', num2str(M(1))];
        dz = 1/M(1);
        z_10 = 0:dz:1;
        FM = ones(1, M(1) - 1);
        tic;
        [F_BFGS_10, A_BFGS_M(1)] = fminunc(@(F) Area_CFD2(F, z_10, F0, F1), FM, options);
        t_BFGS_M(1) = toc;
        F_BFGS_10 = [F0 F_BFGS_10 F1];
        
        % Para 50 puntos:
        M2 = ['N = ', num2str(M(2))];
        dz = 1/M(2);
        z_50 = 0:dz:1;
        FM = ones(1, M(2) - 1);
        tic;
        [F_BFGS_50, A_BFGS_M(2)] = fminunc(@(F) Area_CFD2(F, z_50, F0, F1), FM, options);
        t_BFGS_M(2) = toc;
        F_BFGS_50 = [F0 F_BFGS_50 F1];        

        % Para 100 puntos:
        M3 = ['N = ', num2str(M(3))];
        dz = 1/M(3);
        z_100 = 0:dz:1;
        FM = ones(1, M(3) - 1);
        tic;
        [F_BFGS_100, A_BFGS_M(3)] = fminunc(@(F) Area_CFD2(F, z_100, F0, F1), FM, options);
        t_BFGS_M(3) = toc;
        F_BFGS_100 = [F0 F_BFGS_100 F1];

        % Para 250 puntos:
        M4 = ['N = ', num2str(M(4))];
        dz = 1/M(4);
        z_250 = 0:dz:1;
        FM = ones(1, M(4) - 1);
        tic;
        [F_BFGS_250, A_BFGS_M(4)] = fminunc(@(F) Area_CFD2(F, z_250, F0, F1), FM, options);
        t_BFGS_M(4) = toc;
        F_BFGS_250 = [F0 F_BFGS_250 F1];
        
            
    % Comparación para distinto número de puntos.  
    count = count + 1;
    figure(count)
    plot(F_BFGS_10, z_10, F_BFGS_50, z_50, F_BFGS_100, z_100, F_BFGS_250, z_250);
    axis([0.8 1 0 1])
    title('Influencia del número de puntos (BFGS)')
    legend(M1, M2, M3, M4)
    xlabel('r')
    ylabel('z')
    grid on
    
    
    % Cambio en la distribución de los puntos.
    % Optimización BFGS para distribución sinusoidal de puntos.
    tic;
    [F_BFGS_sin, A_BFGS_sin] = fminunc(@(F) Area_CFD2(F, z_sin, F0, F1), F, options);
    t_BFGS_sin = toc;
    F_BFGS_sin = [F0 F_BFGS_sin F1];
        
        % Ploteado para la distribución sinusoidal.
        count = count + 1;
        figure(count)
        plot(F_BFGS_sin, z_sin);
        axis([0.7 1, 0 1])
        title('Optimización por método BFGS: distribución sinusoidal de puntos')
        xlabel('r')
        ylabel('z')
        grid on
       
        % Error de la distribución sinusoidal.    
        E_BFGS_sin = abs(F_BFGS_sin - FA_sin);
        count = count + 1;
        figure(count)
        plot(E_BFGS_sin, z_sin);
        title('Error del método BFGS: distribución sinusoidal')
        xlabel('Error BFGS')
        ylabel('z')
        grid on
    
    % Comparación de las distintas distribuciones de puntos.
    count = count + 1;
    figure(count)
    plot(F_BFGS_CFD, z, F_BFGS_sin, z_sin);
    axis([0.8 1 0 1])
    title('Influencia de la distribución de los puntos (BFGS)')
    legend('Distribución equiespaciada', 'Distribución sinusoidal')
    xlabel('r')
    ylabel('z')
    grid on
    
    
    % Cambio en la discretización de las derivadas.
    tic;
    [F_BFGS_AtFD, A_BFGS_AtFD] = fminunc(@(F) Area_AtFD3(F, z, F0, F1), F, options);
    t_BFGS_AtFD = toc;

        % Ploteado del resultado obtenido mediante optimización BFGS.
        count = count + 1;
        figure(count)
        F_BFGS_AtFD = [F0 F_BFGS_AtFD F1];
        plot(F_BFGS_AtFD, z);
        axis([0.7 1, 0 1])
        title('Optimización por método BFGS: DFs atrasadas')
        xlabel('r')
        ylabel('z')
        grid on
    
        % Cálculo del error del método BFGS, DFs atrasadas.
        E_BFGS_AtFD = abs(F_BFGS_AtFD - FA);
    
        count = count + 1;
        figure(count)
        plot(E_BFGS_AtFD, z);
        title('Error del método BFGS: DFs atrasadas')
        xlabel('Error BFGS')
        ylabel('z')
        grid on

    % Comparación entre discretizaciones.
    count = count + 1;
    figure(count)
    plot(F_BFGS_CFD, z, F_BFGS_AtFD, z);
    axis([0.9 1, 0 0.1])
    title('Influencia de la discretización (BFGS)')
    legend('DFs Centradas', 'DFs Atrasadas')
    xlabel('r')
    ylabel('z')
    grid on
    
    

%% 4. OPTIMIZACIÓN CON MÉTODO STEEPEST DESCENT

% Método de optimización basado en el gradiente.
% Desarrollado por Peter Debye, en 1909.

% Opciones de optimización para la búsqueda de mínimos mediante Steepest Descent.
    options = optimoptions('fminunc');
    options = optimoptions(options, 'Display', 'off');
    %options = optimoptions(options, 'Display', 'iter');
    options = optimoptions(options, 'MaxFunctionEvaluations', 3e5);           
    options = optimoptions(options, 'MaxIterations', 5e4);               
    options = optimoptions(options, 'FunctionTolerance', 1e-8);                 
    options = optimoptions(options, 'OptimalityTolerance', 1e-8);
    options = optimoptions(options, 'StepTolerance', 1e-8);
    options = optimoptions(options, 'Algorithm', 'quasi-newton');    
    options = optimoptions(options, 'HessUpdate', 'steepdesc');
    
% Optimización Steepest Descent
    tic;
    [F_SD_CFD, A_SD_CFD] = fminunc(@(F) Area_CFD2(F, z, F0, F1), F, options);
    t_SD_CFD = toc;
    
    % Cálculo del volumen Steepest Descent
    F_SD_CFD = [F0 F_SD_CFD F1];
    V_SD_CFD = trapz(z, F_SD_CFD.^2);
    
% Ploteado del resultado obtenido mediante optimización Steepest Descent.
    count = count + 1;
    figure(count)
    plot(F_SD_CFD, z);
    axis([0.7 1, 0 1])
    title('Optimización por método Steepest Descent')
    xlabel('r')
    ylabel('z')
    grid on
    
% Cálculo del error del método Steepest Descent, DFs centradas.
    E_SD_CFD = abs(F_SD_CFD - FA);
    count = count + 1;
    figure(count)
    plot(E_SD_CFD, z);
    title('Error del método Steepest Descent')
    xlabel('Error Steepest Descent')
    ylabel('z')
    grid on
    
    
% Variación de los parámetros del problema.
    
   % Cambio en los radios de los anillos (condiciones de contorno).
   % Vector de valores.
   r0_anillos = [0.8 0.9 1 1.25 1.5];
        
        % Para F0 = 0.8:
        F01 = ['F_0 = ', num2str(r0_anillos(1))];
        F0_CC = r0_anillos(1);
        F1_CC = r0_anillos(1);
        F_CC(1:N-1) = r0_anillos(1);
        tic;
        [F_SD_temp, A_SD_CC(1)] = fminunc(@(F) Area_CFD2(F, z, F0_CC, F1_CC), F_CC, options);
        t_SD_CC(1) = toc;
        F_SD_CC(1,:) = [F0_CC F_SD_temp F1_CC];
           
        % Para F0 = 0.9:
        F02 = ['F_0 = ', num2str(r0_anillos(2))];
        F0_CC = r0_anillos(2);
        F1_CC = r0_anillos(2);
        F_CC(1:N-1) = r0_anillos(2);
        tic;
        [F_SD_temp, A_SD_CC(2)] = fminunc(@(F) Area_CFD2(F, z, F0_CC, F1_CC), F_CC, options); 
        t_SD_CC(2) = toc;
        F_SD_CC(2,:) = [F0_CC F_SD_temp F1_CC];   

        % Para F0 = 1:
        F03 = ['F_0 = ', num2str(r0_anillos(3))];
        F0_CC = r0_anillos(3);
        F1_CC = r0_anillos(3);
        F_CC(1:N-1) = r0_anillos(3);
        tic;
        [F_SD_temp, A_SD_CC(3)] = fminunc(@(F) Area_CFD2(F, z, F0_CC, F1_CC), F_CC, options); 
        t_SD_CC(3) = toc;
        F_SD_CC(3,:) = [F0_CC F_SD_temp F1_CC];
        
        % Para F0 = 1.25:
        F04 = ['F_0 = ', num2str(r0_anillos(4))];
        F0_CC = r0_anillos(4);
        F1_CC = r0_anillos(4);
        F_CC(1:N-1) = r0_anillos(4);
        tic;
        [F_SD_temp, A_SD_CC(4)] = fminunc(@(F) Area_CFD2(F, z, F0_CC, F1_CC), F_CC, options); 
        t_SD_CC(4) = toc;
        F_SD_CC(4,:) = [F0_CC F_SD_temp F1_CC];
        
        % Para F0 = 1.5:
        F05 = ['F_0 = ', num2str(r0_anillos(5))];
        F0_CC = r0_anillos(5);
        F1_CC = r0_anillos(5);
        F_CC(1:N-1) = r0_anillos(5);
        tic;
        [F_SD_temp, A_SD_CC(5)] = fminunc(@(F) Area_CFD2(F, z, F0_CC, F1_CC), F_CC, options); 
        t_SD_CC(5) = toc;
        F_SD_CC(5,:) = [F0_CC F_SD_temp F1_CC];
        
    % Comparación para las distintas condiciones de contorno.
    count = count + 1;
    figure(count)
    plot(F_SD_CC(1,:), z, F_SD_CC(2,:), z, F_SD_CC(3,:), z, F_SD_CC(4,:), z, F_SD_CC(5,:), z);
    axis([0 1.5 0 1])
    title('Influencia de las condiciones de contorno (Steepest Descent)')
    legend(F01, F02, F03, F04, F05)
    xlabel('r')
    ylabel('z')
    grid on
    
    
    % Cambio en el número de puntos.
    M = [10, 50, 100, 250];
    
        % Para 10 puntos:
        M1 = ['N = ', num2str(M(1))];
        dz = 1/M(1);
        z_10 = 0:dz:1;
        FM = ones(1, M(1) - 1);
        tic;
        [F_SD_10, A_SD_M(1)] = fminunc(@(F) Area_CFD2(F, z_10, F0, F1), FM, options);
        t_SD_M(1) = toc;
        F_SD_10 = [F0 F_SD_10 F1];
        
        % Para 50 puntos:
        M2 = ['N = ', num2str(M(2))];
        dz = 1/M(2);
        z_50 = 0:dz:1;
        FM = ones(1, M(2) - 1);
        tic;
        [F_SD_50, A_SD_M(2)] = fminunc(@(F) Area_CFD2(F, z_50, F0, F1), FM, options);
        t_SD_M(2) = toc;
        F_SD_50 = [F0 F_SD_50 F1];        
    
        % Para 100 puntos:
        M3 = ['N = ', num2str(M(3))];
        dz = 1/M(3);
        z_100 = 0:dz:1;
        FM = ones(1, M(3) - 1);
        tic;
        [F_SD_100, A_SD_M(3)] = fminunc(@(F) Area_CFD2(F, z_100, F0, F1), FM, options);
        t_SD_M(3) = toc;
        F_SD_100 = [F0 F_SD_100 F1];

        % Para 250 puntos:
        M4 = ['N = ', num2str(M(4))];
        dz = 1/M(4);
        z_250 = 0:dz:1;
        FM = ones(1, M(4) - 1);
        tic;
        [F_SD_250, A_SD_M(4)] = fminunc(@(F) Area_CFD2(F, z_250, F0, F1), FM, options);
        t_SD_M(4) = toc;
        F_SD_250 = [F0 F_SD_250 F1];
     
        
    % Comparación para distinto número de puntos.
    count = count + 1;
    figure(count)
    plot(F_SD_10, z_10, F_SD_50, z_50, F_SD_100, z_100, F_SD_250, z_250);
    axis([0.8 1 0 1])
    title('Influencia del número de puntos (Steepest Descent)')
    legend(M1, M2, M3, M4)
    xlabel('r')
    ylabel('z')
    grid on

    
    % Cambio en la distribución de los puntos.    
    % Optimización SD para distribución sinusoidal de puntos.
    tic;
    [F_SD_sin, A_SD_sin] = fminunc(@(F) Area_CFD2(F, z_sin, F0, F1), F, options);
    t_SD_sin = toc;
    F_SD_sin = [F0 F_SD_sin F1];
        
        % Ploteado para la distribución sinusoidal.
        count = count + 1;
        figure(count)
        plot(F_SD_sin, z_sin);
        axis([0.7 1, 0 1])
        title('Optimización por método Steepest Descent: distribución sinusoidal de puntos')
        xlabel('r')
        ylabel('z')
        grid on
       
        % Error de la distribución sinusoidal.
        E_SD_sin = abs(F_SD_sin - FA_sin);
        count = count + 1;
        figure(count)
        plot(E_SD_sin, z_sin);
        title('Error del método Steepest Descent: distribución sinusoidal')
        xlabel('Error Steepest Descent')
        ylabel('z')
        grid on
    
    % Comparación de las distintas distribuciones de puntos.
    count = count + 1;
    figure(count)
    plot(F_SD_CFD, z, F_SD_sin, z_sin);
    axis([0.8 1 0 1])
    title('Influencia de la distribución de los puntos (Steepest Descent)')
    legend('Distribución equiespaciada', 'Distribución sinusoidal')
    xlabel('r')
    ylabel('z')
    grid on
    

    % Cambio en la discretización de las derivadas.
    tic;
    [F_SD_AtFD, A_SD_AtFD] = fminunc(@(F) Area_AtFD3(F, z, F0, F1), F, options);
    t_SD_AtFD = toc;

        % Ploteado del resultado obtenido mediante optimización Steepest Descent.
        count = count + 1;
        figure(count)
        F_SD_AtFD = [F0 F_SD_AtFD F1];
        plot(F_SD_AtFD, z);
        axis([0.7 1, 0 1])
        title('Optimización por método Steepest Descent: DFs atrasadas')
        xlabel('r')
        ylabel('z')
        grid on
    
        % Cálculo del error del método Steepest Descent, DFs atrasadas.
        E_SD_AtFD = abs(F_SD_AtFD - FA);
    
        count = count + 1;
        figure(count)
        plot(E_SD_AtFD, z);
        title('Error del método Steepest Descent: DFs atrasadas')
        xlabel('Error Steepest Descent')
        ylabel('z')
        grid on

    % Comparación entre discretizaciones.
    count = count + 1;
    figure(count)
    plot(F_SD_CFD, z, F_SD_AtFD, z);
    axis([0.7 1, 0 1])
    title('Influencia de la discretización (Steepest Descent)')
    legend('DFs Centradas', 'DFs Atrasadas')
    xlabel('r')
    ylabel('z')
    grid on    

    
    
%% 5. OPTIMIZACIÓN CON MÉTODOS HEURÍSTICOS
% ALGORITMO GENÉTICO (GA)
% Desarrollado por John Holland (1975) y mejorado por David E. Goldberg (1989).

% Se reduce el número de puntos a evaluar para mejorar el coste
% computacional y agilizar la obtención de resultados, a costa
% de perder algo de precisión en estos.

    N_GA = 20;
    delta_GA = 1/N_GA;
    z_GA = 0:delta_GA:1;

% Selección de los límites
    UB = ones(1, N_GA-1);           % upper bound
    LB = ones(1, N_GA-1) * 0.75;    % lower bound

% Opciones de optimización para la búsqueda de mínimos mediante GA.
    options = optimoptions('GA');
    options = optimoptions(options, 'Display', 'off');
    %options = optimoptions(options, 'Display', 'iter');

% Optimización GA
    tic
    [F_GA_CFD, A_GA_CFD] = ga(@(F) Area_CFD2(F, z_GA, F0, F1), N_GA - 1, [], [], [], [], LB, UB, [], [], options);
    t_GA_CFD = toc;
    
    % Cálculo del volumen GA
    F_GA_CFD = [F0 F_GA_CFD F1];
    V_GA_CFD = trapz(z_GA, F_GA_CFD.^2);
    
% Ploteado del resultado obtenido mediante optimización GA.
    count = count + 1;
    figure(count)
    plot(F_GA_CFD, z_GA)
    axis([0.7 1, 0 1])
    title('Optimización por método GA')
    xlabel('r')
    ylabel('z')
    grid on
    
% Cálculo del error del método GA. DFs centradas.
    % Obtención de la solución analítica para N_GA puntos.
    [FA_GA, A_FA_GA, V_FA_GA] = Analytical_Function(z_GA, F0, F1);
    
    % Obtención del error propiamente.
    E_GA_CFD = abs(FA_GA - F_GA_CFD);
    count = count + 1;
    figure(count)
    plot(E_GA_CFD, z_GA)
    title('Error del método GA')
    xlabel('Error GA')
    ylabel('z')
    grid on
    
    
% Variación de los parámetros del problema.
    
   % Cambio en los radios de los anillos (condiciones de contorno).
   % Vector de valores.
   r0_anillos_GA = [0.8 0.9 1 1.25 1.5];
        
        % Para F0 = 0.8:
        UB_CC = ones(1, N_GA - 1) * 0.8;
        LB_CC = ones(1, N_GA - 1) * 0.4;
        F01 = ['F_0 = ', num2str(r0_anillos_GA(1))];
        F0_CC = r0_anillos_GA(1);
        F1_CC = r0_anillos_GA(1);
        tic;
        [F_GA_temp, A_GA_CC(1)] = ga(@(F) Area_CFD2(F, z_GA, F0_CC, F1_CC), N_GA - 1, [], [], [], [], LB_CC, UB_CC, [], [], options);
        t_GA_CC(1) = toc;
        F_GA_CC(1,:) = [F0_CC F_GA_temp F1_CC];
           
        % Para F0 = 0.9:
        UB_CC = ones(1, N_GA - 1) * 0.9;
        LB_CC = ones(1, N_GA - 1) * 0.5;
        F02 = ['F_0 = ', num2str(r0_anillos_GA(2))];
        F0_CC = r0_anillos_GA(2);
        F1_CC = r0_anillos_GA(2);
        tic;
        [F_GA_temp, A_GA_CC(2)] = ga(@(F) Area_CFD2(F, z_GA, F0_CC, F1_CC), N_GA - 1, [], [], [], [], LB_CC, UB_CC, [], [], options);
        t_GA_CC(2) = toc;
        F_GA_CC(2,:) = [F0_CC F_GA_temp F1_CC];   

        % Para F0 = 1:
        UB_CC = ones(1, N_GA - 1) * 1;
        LB_CC = ones(1, N_GA - 1) * 0.75;
        F03 = ['F_0 = ', num2str(r0_anillos_GA(3))];
        F0_CC = r0_anillos_GA(3);
        F1_CC = r0_anillos_GA(3);
        tic;
        [F_GA_temp, A_GA_CC(3)] = ga(@(F) Area_CFD2(F, z_GA, F0_CC, F1_CC), N_GA - 1, [], [], [], [], LB_CC, UB_CC, [], [], options);
        t_GA_CC(3) = toc;
        F_GA_CC(3,:) = [F0_CC F_GA_temp F1_CC];
        
        % Para F0 = 1.25:
        UB_CC = ones(1, N_GA - 1) * 1.25;
        LB_CC = ones(1, N_GA - 1) * 0.8;
        F04 = ['F_0 = ', num2str(r0_anillos_GA(4))];
        F0_CC = r0_anillos_GA(4);
        F1_CC = r0_anillos_GA(4);
        tic;
        [F_GA_temp, A_GA_CC(4)] = ga(@(F) Area_CFD2(F, z_GA, F0_CC, F1_CC), N_GA - 1, [], [], [], [], LB_CC, UB_CC, [], [], options);
        t_GA_CC(4) = toc;
        F_GA_CC(4,:) = [F0_CC F_GA_temp F1_CC];
        
        % Para F0 = 1.5:
        UB_CC = ones(1, N_GA - 1) * 1.5;
        LB_CC = ones(1, N_GA - 1) * 1;
        F05 = ['F_0 = ', num2str(r0_anillos_GA(5))];
        F0_CC = r0_anillos_GA(5);
        F1_CC = r0_anillos_GA(5);
        tic;
        [F_GA_temp, A_GA_CC(5)] = ga(@(F) Area_CFD2(F, z_GA, F0_CC, F1_CC), N_GA - 1, [], [], [], [], LB_CC, UB_CC, [], [], options);
        t_GA_CC(5) = toc;
        F_GA_CC(5,:) = [F0_CC F_GA_temp F1_CC];
        
    % Comparación para las distintas condiciones de contorno.
    count = count + 1;
    figure(count)
    plot(F_GA_CC(1,:), z_GA, F_GA_CC(2,:), z_GA, F_GA_CC(3,:), z_GA, F_GA_CC(4,:), z_GA, F_GA_CC(5,:), z_GA);
    axis([0 1.5 0 1])
    title('Influencia de las condiciones de contorno (GA)')
    legend(F01, F02, F03, F04, F05)
    xlabel('r')
    ylabel('z')
    grid on    

    
    % Cambio en el número de puntos.
    m = [5, 15, 30, 60];
    
        % Para 5 puntos:
        LBm = ones(1, m(1) - 1) * 0.75;
        UBm = ones(1, m(1) - 1);
        m1 = ['N = ', num2str(m(1))];
        dz = 1/m(1);
        z_GA_5 = 0:dz:1;
        tic;
        [F_GA_5, A_GA_m(1)] = ga(@(F) Area_CFD2(F, z_GA_5, F0, F1), m(1) - 1, [], [], [], [], LBm, UBm, [], [], options);
        t_GA_m(1) = toc;
        F_GA_5 = [F0 F_GA_5 F1];
        
        % Para 15 puntos:
        LBm = ones(1, m(2) - 1) * 0.75;
        UBm = ones(1, m(2) - 1);
        m2 = ['N = ', num2str(m(2))];
        dz = 1/m(2);
        z_GA_15 = 0:dz:1;
        tic;
        [F_GA_15, A_GA_m(2)] = ga(@(F) Area_CFD2(F, z_GA_15, F0, F1), m(2) - 1, [], [], [], [], LBm, UBm, [], [], options);
        t_GA_m(2) = toc;
        F_GA_15 = [F0 F_GA_15 F1];        

        % Para 30 puntos:
        LBm = ones(1, m(3) - 1) * 0.75;
        UBm = ones(1, m(3) - 1);
        m3 = ['N = ', num2str(m(3))];
        dz = 1/m(3);
        z_GA_30 = 0:dz:1;
        tic;
        [F_GA_30, A_GA_m(3)] = ga(@(F) Area_CFD2(F, z_GA_30, F0, F1), m(3) - 1, [], [], [], [], LBm, UBm, [], [], options);
        t_GA_m(3) = toc;
        F_GA_30 = [F0 F_GA_30 F1];

        % Para 60 puntos:
        LBm = ones(1, m(4) - 1) * 0.75;
        UBm = ones(1, m(4) - 1);
        m4 = ['N = ', num2str(m(4))];
        dz = 1/m(4);
        z_GA_60 = 0:dz:1;
        tic;
        [F_GA_60, A_GA_m(4)] = ga(@(F) Area_CFD2(F, z_GA_60, F0, F1), m(4) - 1, [], [], [], [], LBm, UBm, [], [], options);
        t_GA_m(4) = toc;
        F_GA_60 = [F0 F_GA_60 F1];
        
    % Comparación para distinto número de puntos.  
    count = count + 1;
    figure(count)
    plot(F_GA_5, z_GA_5, F_GA_15, z_GA_15, F_GA_30, z_GA_30, F_GA_60, z_GA_60);
    axis([0.8 1 0 1])
    title('Influencia del número de puntos (GA)')
    legend(m1, m2, m3, m4)
    xlabel('r')
    ylabel('z')
    grid on
    
    
    % Cambio en la distribución de los puntos.
    % Distribución sinusoidal simple.
    z_GA_sin = sin((0:delta_GA:1) * pi / 2);
    
    % Optimización GA para distribución sinusoidal de puntos.
    tic;
    [F_GA_sin, A_GA_sin] = ga(@(F) Area_CFD2(F, z_GA_sin, F0, F1), N_GA - 1, [], [], [], [], LB, UB, [], [], options);
    t_GA_sin = toc;
    F_GA_sin = [F0 F_GA_sin F1];
        
        % Ploteado para la distribución sinusoidal.
        count = count + 1;
        figure(count)
        plot(F_GA_sin, z_GA_sin);
        axis([0.7 1, 0 1])
        title('Optimización por método GA: distribución sinusoidal de puntos')
        xlabel('r')
        ylabel('z')
        grid on
       
        % Error de la distribución sinusoidal.
        % Obtención de la solución analítica para la distribución sinusoidal.
        [FA_GA_sin, A_FA_GA_sin, V_FA_GA_sin] = Analytical_Function(z_GA_sin, F0, F1);
        
        E_GA_sin = abs(F_GA_sin - FA_GA_sin);
        count = count + 1;
        figure(count)
        plot(E_GA_sin, z_GA_sin);
        title('Error del método GA: distribución sinusoidal')
        xlabel('Error GA')
        ylabel('z')
        grid on
    
    % Comparación de las distintas distribuciones de puntos.
    count = count + 1;
    figure(count)
    plot(F_GA_CFD, z_GA, F_GA_sin, z_GA_sin);
    axis([0.8 1 0 1])
    title('Influencia de la distribución de los puntos (GA)')
    legend('Distribución equiespaciada', 'Distribución sinusoidal')
    xlabel('r')
    ylabel('z')
    grid on

    
    % Cambio en la discretización de las derivadas.
    tic;
    [F_GA_AtFD, A_GA_AtFD] = ga(@(F) Area_AtFD3(F, z_GA, F0, F1), N_GA - 1, [], [], [], [], LB, UB, [], [], options);
    t_GA_AtFD = toc;

        % Ploteado del resultado obtenido mediante optimización GA.
        count = count + 1;
        figure(count)
        F_GA_AtFD = [F0 F_GA_AtFD F1];
        plot(F_GA_AtFD, z_GA);
        axis([0.7 1, 0 1])
        title('Optimización por método GA: DFs atrasadas')
        xlabel('r')
        ylabel('z')
        grid on
    
        % Cálculo del error del método GA, DFs atrasadas.
        E_GA_AtFD = abs(F_GA_AtFD - FA_GA);
    
        count = count + 1;
        figure(count)
        plot(E_GA_AtFD, z_GA);
        title('Error del método GA: DFs atrasadas')
        xlabel('Error GA')
        ylabel('z')
        grid on

    % Comparación entre discretizaciones.
    count = count + 1;
    figure(count)
    plot(F_GA_CFD, z_GA, F_GA_AtFD, z_GA);
    axis([0.7 1, 0 1])
    title('Influencia de la discretización (GA)')
    legend('DFs Centradas', 'DFs Atrasadas')
    xlabel('r')
    ylabel('z')
    grid on
    
    
    
%% 6. ECUACIÓN DE EULER
% Apartado opcional.
% Resuelve numéricamente la ecuación: F * F'' = 1 + (F')^2.
    
    euler_options = optimoptions('fsolve');
    euler_options = optimoptions(euler_options, 'Display', 'off');
    %euler_options = optimoptions(euler_options, 'Display', 'iter');
    euler_options = optimoptions(euler_options, 'Algorithm', 'levenberg-marquardt');
    % Se ha comprobado que el problema no admite el algoritmo 
    % predeterminado 'trust-region-dogleg'.
    
% Obtención de la solución por la Ecuación de Euler.
    tic;
    F_Euler = fsolve(@(F) Euler_Equation(F, z, F0, F1), F, euler_options);
    t_Euler = toc;
    F_Euler = [F0 F_Euler F1];
    
% Cálculo del área y volumen.
    dF_Euler = zeros(1, N+1);
    % Derivada de la función.
    for i = 1:N+1
        % extremo inferior obliga a utilizar derivada adelantada.
        if i == 1
            dF_Euler(i) = (F_Euler(i+1) - F_Euler(i)) / (z(i+1) - z(i));
        % extremo superior obliga a utilizar derivada atrasada.
        elseif i == N+1
            dF_Euler(i) = (F_Euler(i) - F_Euler(i-1)) / (z(i) - z(i-1));
        else    
            dF_Euler(i) = (F_Euler(i+1) - F_Euler(i-1)) / (z(i+1) - z(i-1));
        end
    end
    
    I_A_Euler = F_Euler .* sqrt(1 + dF_Euler.^2);
    I_V_Euler = F_Euler.^2;
    A_Euler = trapz(z, I_A_Euler);
    V_Euler = trapz(z, I_V_Euler);
    
    
% Ploteado de la solución obtenida por este método.
    count = count + 1;
    figure(count)
    plot(F_Euler, z)
    axis([0.7 1, 0 1])
    title('Solución por la Ecuación de Euler')
    xlabel('r')
    ylabel('z')
    grid on

% Cálculo y ploteado del error asociado a la solución por Euler.
    E_Euler = abs(F_Euler - FA);
    
    count = count + 1;
    figure(count)
    plot(E_Euler, z)
    title('Error de la Ecuación de Euler')
    xlabel('Error Euler')
    ylabel('z')
    grid on
    
    
%% 7. PROBLEMA CON RESTRICCIÓN DE VOLUMEN

% Métodos basados en el gradiente, con restricciones.
    
    % Opciones para la optimización con restricciones.
        c_options = optimoptions('fmincon');
        c_options = optimoptions(c_options, 'Display', 'off');
        %c_options = optimoptions(c_options, 'Display', 'iter');
        c_options = optimoptions(c_options, 'Algorithm', 'SQP');
        c_options = optimoptions(c_options, 'HessianApproximation', 'BFGS');
    % SQP (Sequential Quadratic Programming) es un método iterativo para
    % optimización no lineal con restricciones.

    % Restricción implementada mediante el comando 'fmincon'.
        % VC = 100% V_FA
        VC = 1;
        VC1 = ('V = V_0');
        VC = VC * V_FA;
        tic;
        [F_SQP_temp, A_SQP(1)] = fmincon(@(F) Area_CFD2(F, z, F0, F1), F, [], [], [], [], [], [], @(F) Volume_Constraint(F, z, F0, F1, VC), c_options);
        t_SQP(1) = toc;
        F_SQP(1, :) = [F0 F_SQP_temp F1];

        % VC = 75% V_FA
        VC = 0.75;
        VC2 = ['V = ', num2str(VC), 'V_0'];
        VC = VC * V_FA;
        tic;
        [F_SQP_temp, A_SQP(2)] = fmincon(@(F) Area_CFD2(F, z, F0, F1), F, [], [], [], [], [], [], @(F) Volume_Constraint(F, z, F0, F1, VC), c_options);
        t_SQP(2) = toc;
        F_SQP(2, :) = [F0 F_SQP_temp F1]; 
    
        % VC = 125% V_FA
        VC = 1.25;
        VC3 = ['V = ', num2str(VC), 'V_0'];
        VC = VC * V_FA;
        tic;
        [F_SQP_temp, A_SQP(3)] = fmincon(@(F) Area_CFD2(F, z, F0, F1), F, [], [], [], [], [], [], @(F) Volume_Constraint(F, z, F0, F1, VC), c_options);
        t_SQP(3) = toc;
        F_SQP(3, :) = [F0 F_SQP_temp F1]; 
    
        % VC = 150% V_FA
        VC = 1.5;
        VC4 = ['V = ', num2str(VC), 'V_0'];
        VC = VC * V_FA;
        tic;
        [F_SQP_temp, A_SQP(4)] = fmincon(@(F) Area_CFD2(F, z, F0, F1), F, [], [], [], [], [], [], @(F) Volume_Constraint(F, z, F0, F1, VC), c_options);
        t_SQP(4) = toc;
        F_SQP(4, :) = [F0 F_SQP_temp F1];  
    
    % Ploteado de las soluciones para las distintas restricciones de
    % volumen. Mediante Sequential Quadratic Programming.
        count = count + 1;
        figure(count)
        plot(F_SQP(1,:), z, F_SQP(2,:), z, F_SQP(3,:), z, F_SQP(4,:), z)
        axis([0.5 1.5, 0 1])
        title('Optimización bajo restricciones: método SQP')
        legend(VC1, VC2, VC3, VC4)
        xlabel('r')
        ylabel('z')
        grid on

% Métodos heurísticos, con restricciones.

    % Métodos muy computacionalmente exigentes. Más vale un número de
    % puntos bajo.
    N_GA_C = 20;
    delta_GA_C = 1/N_GA_C;
    z_GA_C = 0:delta_GA_C:1;
    
    % Selección de los límites.
    LB = ones(1, N_GA_C - 1) * 0;
    UB = ones(1, N_GA_C - 1) * 1.5;
    
    c_options = optimoptions('GA');
    c_options = optimoptions(c_options, 'Display', 'off');
    %c_options = optimoptions(c_options, 'Display', 'iter');
    
    % El comando 'ga' ('genetic algorithm') da la posibilidad de introducir
    % las ligaduras deseadas.
        % VC = 100% V_FA
        VC = 1;
        VC1 = ('V = V_0');
        VC = VC * V_FA;
        tic;
        [F_GA_temp, A_GA_C(1)] = ga(@(F) Area_CFD2(F, z_GA_C, F0, F1), N_GA_C - 1, [], [], [], [], LB, UB, @(F) Volume_Constraint(F, z_GA_C, F0, F1, VC), [], options);
        t_GA_C(1) = toc;
        F_GA_C(1,:) = [F0 F_GA_temp F1];
        
        % VC = 75% V_FA
        VC = 0.75;
        VC2 = ['V = ', num2str(VC), 'V_0'];
        VC = VC * V_FA;
        tic;
        [F_GA_temp, A_GA_C(2)] = ga(@(F) Area_CFD2(F, z_GA_C, F0, F1), N_GA_C - 1, [], [], [], [], LB, UB, @(F) Volume_Constraint(F, z_GA_C, F0, F1, VC), [], options);
        t_GA_C(2) = toc;
        F_GA_C(2,:) = [F0 F_GA_temp F1];
        
        % VC = 125% V_FA
        VC = 1.25;
        VC3 = ['V = ', num2str(VC), 'V_0'];
        VC = VC * V_FA;
        tic;
        [F_GA_temp, A_GA_C(3)] = ga(@(F) Area_CFD2(F, z_GA_C, F0, F1), N_GA_C - 1, [], [], [], [], LB, UB, @(F) Volume_Constraint(F, z_GA_C, F0, F1, VC), [], options);
        t_GA_C(3) = toc;
        F_GA_C(3,:) = [F0 F_GA_temp F1];
        
        % VC = 150% V_FA
        VC = 1.5;
        VC4 = ['V = ', num2str(VC), 'V_0'];
        VC = VC * V_FA;
        tic;
        [F_GA_temp, A_GA_C(4)] = ga(@(F) Area_CFD2(F, z_GA_C, F0, F1), N_GA_C - 1, [], [], [], [], LB, UB, @(F) Volume_Constraint(F, z_GA_C, F0, F1, VC), [], options);
        t_GA_C(4) = toc;
        F_GA_C(4,:) = [F0 F_GA_temp F1];
        
    % Ploteado de las soluciones para las distintas restricciones de
    % volumen. Mediante Genetic Algorithm.
    
        count = count + 1;
        figure(count)
        plot(F_GA_C(1,:), z_GA_C, F_GA_C(2,:), z_GA_C, F_GA_C(3,:), z_GA_C, F_GA_C(4,:), z_GA_C)
        axis([0.5 1.5, 0 1])
        title('Optimización bajo restricciones: método GA')
        legend(VC1, VC2, VC3, VC4)
        xlabel('r')
        ylabel('z')
        grid on
        
        
%% 8. PROBLEMA CON LOS RADIOS DESIGUALES

% Opciones de optimización para la búsqueda de mínimos mediante BFGS.
    options = optimoptions('fminunc');
    options = optimoptions(options, 'Display', 'off');
    %options = optimoptions(options, 'Display', 'iter');
    options = optimoptions(options, 'MaxFunctionEvaluations', 2e6);
    options = optimoptions(options, 'MaxIterations', 1e4);
    options = optimoptions(options, 'FunctionTolerance', 1e-8);
    options = optimoptions(options, 'OptimalityTolerance', 1e-8);
    options = optimoptions(options, 'StepTolerance', 1e-8);
    options = optimoptions(options, 'Algorithm', 'quasi-newton');
    options = optimoptions(options, 'ObjectiveLimit', 1e-6);
    options = optimoptions(options, 'HessUpdate', 'BFGS');
    
    
% Se estudia el problema para distintas perturbaciones.
EPS = [0.1 0.3 0.5];        % desigualdad en los radios de los anillos

    % EPS = 0.0
    % Hecho en línea 320
    EPS0 = ('eps = 0');
    t_BFGS_EPS(1) = t_BFGS_CFD;
    
    % EPS = 0.1
    EPS1 = ['eps = ', num2str(EPS(1))];
    F0_EPS = r0 - EPS(1);
    F1_EPS = r0 + EPS(1);    
    F_EPS(1:N-1) = F0_EPS + z(1:N-1) * (F1_EPS - F0_EPS); 
    
    tic;
    [F_BFGS_temp, A_BFGS_EPS(1)] = fminunc(@(F) Area_CFD2(F, z, F0_EPS, F1_EPS), F_EPS, options);
    t_BFGS_EPS(2) = toc;
    F_BFGS_EPS(1,:) = [F0_EPS F_BFGS_temp F1_EPS];
   
    % EPS = 0.3
    EPS2 = ['eps = ', num2str(EPS(2))];
    F0_EPS = r0 - EPS(2);
    F1_EPS = r0 + EPS(2);    
    F_EPS(1:N-1) = F0_EPS + z(1:N-1) * (F1_EPS - F0_EPS); 
    
    tic;
    [F_BFGS_temp, A_BFGS_EPS(2)] = fminunc(@(F) Area_CFD2(F, z, F0_EPS, F1_EPS), F_EPS, options);
    t_BFGS_EPS(3) = toc;
    F_BFGS_EPS(2,:) = [F0_EPS F_BFGS_temp F1_EPS];
    
    % EPS = 0.5
    EPS3 = ['eps = ', num2str(EPS(3))];
    F0_EPS = r0 - EPS(3);
    F1_EPS = r0 + EPS(3);    
    F_EPS(1:N-1) = F0_EPS + z(1:N-1) * (F1_EPS - F0_EPS); 
    
    tic;
    [F_BFGS_temp, A_BFGS_EPS(3)] = fminunc(@(F) Area_CFD2(F, z, F0_EPS, F1_EPS), F_EPS, options);
    t_BFGS_EPS(4) = toc;
    F_BFGS_EPS(3,:) = [F0_EPS F_BFGS_temp F1_EPS];
    
% Ploteado de las soluciones para distintas perturbaciones en los radios de
% los anillos. Método BFGS, basado en el gradiente.

    count = count + 1;
    figure(count)
    plot(F_BFGS_CFD, z, F_BFGS_EPS(1,:), z, F_BFGS_EPS(2,:), z, F_BFGS_EPS(3,:), z)
    axis([0 2, 0 1])
    title('Optimización para radios distintos: método BFGS')
    legend(EPS0, EPS1, EPS2, EPS3)
    xlabel('r')
    ylabel('z')
    grid on
   
    
% Métodos heurísticos, para discos desiguales.

% Métodos muy computacionalmente exigentes. Más vale un número de
% puntos bajo.
N_GA_EPS = 20;
delta_GA_EPS = 1/N_GA_EPS;
z_GA_EPS = 0:delta_GA_EPS:1;
    
% Selección de los límites.
LB = ones(1, N_GA_EPS - 1) * 0.4;
UB = ones(1, N_GA_EPS - 1) * 1.6;

% Optimización mediante método GA.
    options = optimoptions('GA');
    options = optimoptions(options, 'Display', 'off');
    %options = optimoptions(options, 'Display', 'iter');
        
% Se estudia el problema para distintas perturbaciones.
EPS = [0.1 0.3 0.5];        % desigualdad en los radios de los anillos

% Inicialización de la matriz de funciones.
F_GA_EPS = zeros(max(size(EPS)), N_GA_EPS+1);

    % EPS = 0.0
    % Hecho en línea 823
    EPS0 = ('eps = 0');
    t_GA_EPS(1) = t_GA_CFD;
    
    % EPS = 0.1
    EPS1 = ['eps = ', num2str(EPS(1))];
    F0_EPS = r0 - EPS(1);
    F1_EPS = r0 + EPS(1);    
    F_EPS(1:N_GA_EPS-1) = F0_EPS + z_GA_EPS(1:N_GA_EPS-1) * (F1_EPS - F0_EPS);
    
    tic;
    [F_GA_temp, A_GA_EPS(1)] = ga(@(F) Area_CFD2(F, z_GA_EPS, F0_EPS, F1_EPS), N_GA_EPS - 1, [], [], [], [], LB, UB, [], [], options);
    t_GA_EPS(2) = toc;
    F_GA_EPS(1,:) = [F0_EPS F_GA_temp F1_EPS];
   
    % EPS = 0.3
    EPS2 = ['eps = ', num2str(EPS(2))];
    F0_EPS = r0 - EPS(2);
    F1_EPS = r0 + EPS(2);    
    F_EPS(1:N_GA_EPS-1) = F0_EPS + z_GA_EPS(1:N_GA_EPS-1) * (F1_EPS - F0_EPS); 
    
    tic;
    [F_GA_temp, A_GA_EPS(2)] = ga(@(F) Area_CFD2(F, z_GA_EPS, F0_EPS, F1_EPS), N_GA_EPS - 1, [], [], [], [], LB, UB, [], [], options);
    t_GA_EPS(3) = toc;
    F_GA_EPS(2,:) = [F0_EPS F_GA_temp F1_EPS];
    
    % EPS = 0.5
    EPS3 = ['eps = ', num2str(EPS(3))];
    F0_EPS = r0 - EPS(3);
    F1_EPS = r0 + EPS(3);    
    F_EPS(1:N_GA_EPS-1) = F0_EPS + z_GA_EPS(1:N_GA_EPS-1) * (F1_EPS - F0_EPS); 
    
    tic;
    [F_GA_temp, A_GA_EPS(3)] = ga(@(F) Area_CFD2(F, z_GA_EPS, F0_EPS, F1_EPS), N_GA_EPS - 1, [], [], [], [], LB, UB, [], [], options);
    t_GA_EPS(4) = toc;
    F_GA_EPS(3,:) = [F0_EPS F_GA_temp F1_EPS];
    
% Ploteado de las soluciones para distintas perturbaciones en los radios de
% los anillos. Método GA.

    count = count + 1;
    figure(count)
    plot(F_GA_CFD, z_GA, F_GA_EPS(1,:), z_GA_EPS, F_GA_EPS(2,:), z_GA_EPS, F_GA_EPS(3,:), z_GA_EPS)
    axis([0 2, 0 1])
    title('Optimización para radios distintos: método GA')
    legend(EPS0, EPS1, EPS2, EPS3)
    xlabel('r')
    ylabel('z')
    grid on
    
    
    
%% FIN DEL CÓDIGO: NO QUEDA YA NADA POR OPTIMIZAR EN ESTA VIDA.
Tiempo_total_de_ejecucion = toc(start);