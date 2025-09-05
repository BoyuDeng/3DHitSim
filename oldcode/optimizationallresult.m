function [optimized_coeffs, optimized_W, totalEnergy, fval, all_solutions] = optimizationallresult(t, uField, vField, wField, dt, p, U, Forcing, customStartPoints)

    % Set bounds for the optimization problem
    lb = [-1 * ones(8, 1); -1 * ones(8, 1); -1 * ones(10, 1); 0.001];  % Lower bounds for coefficients and custom lower bound for W
    ub = [1 * ones(8, 1); 1 * ones(8, 1); 1 * ones(10, 1); 10];  % Upper bounds for coefficients and custom upper bound for W

    % Define optimization options using the interior-point algorithm
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');

    % Define the optimization problem for fmincon
    problem = createOptimProblem('fmincon', 'objective', ...
        @(vars) COT14(vars(1:26), t, vars(27), uField, vField, wField, dt, p, U, Forcing), ...
        'x0', zeros(27, 1), 'lb', lb, 'ub', ub, 'options', options);

<<<<<<< HEAD
    % Validate the custom starting points
    if size(customStartPoints, 2) ~= 27
        error('customStartPoints must have 27 columns (26 coefficients + 1 W).');
    end
=======
    % Create a GlobalSearch object to perform the optimization
    gs = GlobalSearch('NumTrialPoints', 500);
>>>>>>> 1d1a4fd6e081401d697792f152ca234f12bad312

    % Create a CustomStartPointSet object
    startPointSet = CustomStartPointSet(customStartPoints);

    % Create a MultiStart object
    ms = MultiStart('UseParallel', true, 'Display', 'off');

    % Run the multi-start optimization
    [result, fval, exitflag, output, solutions] = run(ms, problem, startPointSet);

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
