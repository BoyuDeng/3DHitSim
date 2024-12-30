function [optimized_coeffs, optimized_W, totalEnergy, fval] = optimization27(t, W, uField, vField, wField, dt, p, U, Forcing)

    % Combine the coefficients and W into a single vector for optimization
    initial_coeffs = zeros(27,1);  % 14 coefficients + 1 for W
    initial_coeffs(27) = W;  % Set the initial guess for W

% Set bounds for the optimization problem
    lb = [-0.1 * ones(8, 1);-1 * ones(8, 1);-0.1 * ones(10, 1); 0.001];  % Lower bounds for coefficients and custom lower bound for W
    ub = [0.1* ones(8, 1);1 * ones(8, 1);0.1 * ones(10, 1); 1.2];    % Upper bounds for coefficients and custom upper bound for W

    % Define optimization options using the interior-point algorithm
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter');

    % Define the optimization problem for fmincon
    problem = createOptimProblem('fmincon', 'objective', ...
        @(vars) COT14(vars(1:26), t, vars(27), uField, vField, wField, dt, p, U, Forcing), ...
        'x0', initial_coeffs, 'lb', lb, 'ub', ub, 'options', options);

    % Create a GlobalSearch object to perform the optimization
    gs = GlobalSearch;

    % Run the global optimization
    [result, fval] = run(gs, problem);

    % Extract the optimized coefficients and W
    optimized_coeffs = result(1:26);
    optimized_W = result(27);

    % Calculate the total energy using the optimized coefficients and W
    totalEnergy = COT14(optimized_coeffs, t, optimized_W, uField, vField, wField, dt, p, U, Forcing);

    % Display the result
    disp('Optimized total energy:');
    disp(totalEnergy);
end
