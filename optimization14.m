function [optimized_coeffs, totalEnergy, fval] = optimization14(t, W, uField, vField, wField, dt, p, U, Forcing)

    % Initial coefficients from A at min_index
    initial_coeffs = zeros(14,1);
    
    % Set bounds for the optimization problem
    lb = -0.1 * ones(12, 1);  % Lower bounds
    ub = 0.1 * ones(12, 1);   % Upper bounds
    
    % Define optimization options using the interior-point algorithm
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter');
    
    % Define the optimization problem for fmincon
    problem = createOptimProblem('fmincon', 'objective', ...
        @(coeffs) COT14(coeffs, t, W, uField, vField, wField, dt, p, U, Forcing), ...
        'x0', initial_coeffs, 'lb', lb, 'ub', ub, 'options', options);
    
    % Create a GlobalSearch object to perform the optimization
    gs = GlobalSearch;
    
    % Run the global optimization
    [result, fval] = run(gs, problem);
    
    % Extract the optimized coefficients
    optimized_coeffs = result;
    
    % Calculate the total energy using the optimized coefficients
    totalEnergy = COT14(optimized_coeffs, t, W, uField, vField, wField, dt, p, U, Forcing);
    
    % Display the result
    disp('Optimized total energy:');
    disp(totalEnergy);
end
