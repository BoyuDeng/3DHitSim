function [optimized_coeffs, optimized_W, totalEnergy, fval] = optimization27(t, W, uField, vField, wField, dt, p, U, Forcing)
    mode = 12;



    % Combine coefficients and W into a single vector
    initial_coeffs = zeros(2*mode+mode/2+3, 1);
    initial_coeffs(end) = W;

    
    % Define bounds
    lb = [-1 * ones(mode, 1); -1 * ones(mode/2, 1); -1 * ones(mode + 2, 1); 5];
    ub = [1 * ones(mode, 1); 1 * ones(mode/2, 1); 1 * ones(mode + 2, 1); 10];



    % Define optimization options using the interior-point algorithm
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'none');

    % Define the optimization problem for fmincon
    problem = createOptimProblem('fmincon', 'objective', ...
        @(vars) COT14(vars(1:end-1), t, vars(end), uField, vField, wField, dt, p, U, Forcing), ...
        'x0', initial_coeffs, 'lb', lb, 'ub', ub, 'options', options);

    % Create a GlobalSearch object to perform the optimization
      gs = GlobalSearch('NumTrialPoints', 10000);

    % Run the global optimization
    [result, fval] = run(gs, problem);

    % Extract the optimized coefficients and W
    optimized_coeffs = result(1:end-1);
    optimized_W = result(end);

    % Calculate the total energy using the optimized coefficients and W
    totalEnergy = COT14(optimized_coeffs, t, optimized_W, uField, vField, wField, dt, p, U, Forcing);

    % Display the result
    % disp('Optimized total energy:');
    % disp(totalEnergy);
end

