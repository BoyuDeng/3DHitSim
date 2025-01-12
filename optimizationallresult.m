function [optimized_coeffs, optimized_W, totalEnergy, fval, all_solutions] = optimizationallresult(t, W, uField, vField, wField, dt, p, U, Forcing)

    % Combine the coefficients and W into a single vector for optimization
    initial_coeffs = zeros(27, 1);  % 14 coefficients + 1 for W
    initial_coeffs(27) = W;  % Set the initial guess for W

    % Set bounds for the optimization problem
    lb = [-1 * ones(8, 1); -10 * ones(8, 1); -1 * ones(10, 1); 0.001];  % Lower bounds for coefficients and custom lower bound for W
    ub = [1 * ones(8, 1); 10 * ones(8, 1); 1 * ones(10, 1); 10];  % Upper bounds for coefficients and custom upper bound for W

    % Define optimization options using the interior-point algorithm
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');

    % Define the optimization problem for fmincon
    problem = createOptimProblem('fmincon', 'objective', ...
        @(vars) COT14(vars(1:26), t, vars(27), uField, vField, wField, dt, p, U, Forcing), ...
        'x0', initial_coeffs, 'lb', lb, 'ub', ub, 'options', options);

    % Create a GlobalSearch object to perform the optimization
    gs = GlobalSearch('NumTrialPoints', 3000);

    % Run the global optimization
    [result, fval, exitflag, output, solutions] = run(gs, problem);

    % Extract the optimized coefficients and W from the best solution
    optimized_coeffs = result(1:26);
    optimized_W = result(27);

    % Calculate the total energy using the optimized coefficients and W
    totalEnergy = COT14(optimized_coeffs, t, optimized_W, uField, vField, wField, dt, p, U, Forcing);

    % Save all solutions (coefficients, W, and fval) in an array
    num_solutions = length(solutions);
    all_solutions = zeros(num_solutions, 28); % 26 coefficients + W + fval
    for i = 1:num_solutions
        all_solutions(i, 1:26) = solutions(i).X(1:26); % Coefficients
        all_solutions(i, 27) = solutions(i).X(27);    % W
        all_solutions(i, 28) = solutions(i).Fval;     % Objective value
    end
end
