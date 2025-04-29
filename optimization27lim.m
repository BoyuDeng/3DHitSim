function [optimized_coeffs, optimized_W, totalEnergy, fval] = optimization27lim(t, W, uField, vField, wField, dt, p, U, Forcing, trialpoints, mode, tau_p, wlow,whigh, lim)
    % mode = mode;
    

    % Combine coefficients and W into a single vector
    initial_coeffs = zeros(2*mode+mode/2+3, 1);
    initial_coeffs(end) = W;

    
    % Define bounds
    lb = [-lim * ones(mode, 1); -lim * ones(mode/2, 1); -lim * ones(mode + 2, 1); wlow];
    ub = [lim * ones(mode, 1); lim * ones(mode/2, 1); lim * ones(mode + 2, 1); whigh];



    % Define optimization options using the interior-point algorithm
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'none');

    % Define the optimization problem for fmincon
    problem = createOptimProblem('fmincon', 'objective', ...
        @(vars) COT14(vars(1:end-1), t, vars(end), uField, vField, wField, dt, p, U, Forcing, tau_p), ...
        'x0', initial_coeffs, 'lb', lb, 'ub', ub, 'options', options);

    % Create a GlobalSearch object to perform the optimization
      gs = GlobalSearch('NumTrialPoints', trialpoints);

    % Run the global optimization
    [result, fval] = run(gs, problem);

    % Extract the optimized coefficients and W
    optimized_coeffs = result(1:end-1);
    optimized_W = result(end);

    % Calculate the total energy using the optimized coefficients and W
    totalEnergy = COT14(optimized_coeffs, t, optimized_W, uField, vField, wField, dt, p, U, Forcing, tau_p);

    % Display the result
    % disp('Optimized total energy:');
    % disp(totalEnergy);
end

