function [optimized_coeffs, optimized_W, totalEnergy, fval] = optMultiStart(t, W, uField, vField, wField, dt, p, U, Forcing)

    % Combine the coefficients and W into a single vector for optimization
    initial_coeffs = zeros(27, 1);  % 26 coefficients + 1 for W
    initial_coeffs(27) = W;         % Set the initial guess for W

    % Set bounds for the optimization problem
    lb = [-1 * ones(8, 1); -1 * ones(8, 1); -1 * ones(10, 1); 0.001];  % Lower bounds
    ub = [1 * ones(8, 1); 1 * ones(8, 1); 1 * ones(10, 1); 1.2];       % Upper bounds

    % Define optimization options for fmincon
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter');

    % Define the optimization problem for fmincon
    problem = createOptimProblem('fmincon', 'objective', ...
        @(vars) COT14(vars(1:26), t, vars(27), uField, vField, wField, dt, p, U, Forcing), ...
        'x0', initial_coeffs, 'lb', lb, 'ub', ub, 'options', options);

    % Define the MultiStart object
    ms = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun', 'all');

    % Generate random starting points for MultiStart
    start_points = RandomStartPointSet('NumStartPoints', 1000);

    % Run the MultiStart optimization
    [result, fval, exitflag, output, solutions] = run(ms, problem, start_points);

    % Extract the optimized coefficients and W
    optimized_coeffs = result(1:26);
    optimized_W = result(27);

    % Calculate the total energy using the optimized coefficients and W
    totalEnergy = COT14(optimized_coeffs, t, optimized_W, uField, vField, wField, dt, p, U, Forcing);

    % Display the result
    disp('Optimized total energy:');
    disp(totalEnergy);
end
