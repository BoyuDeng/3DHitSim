function [optimized_coeffs, optimized_W, totalEnergy, fval, initial_coeffs] = optimization_annealing(t, W, uField, vField, wField, dt, p, U, Forcing)
    
    mode = 8;



    % Combine coefficients and W into a single vector
    initial_coeffs = randn(2*mode + mode/2 + 3, 1);
    initial_coeffs(end) = W;  % Keep W as given

   


    % Define bounds
    lb = [-1 * ones(mode, 1); -1 * ones(mode/2, 1); -1 * ones(mode + 2, 1); 1];
    ub = [1 * ones(mode, 1); 1 * ones(mode/2, 1); 1 * ones(mode + 2, 1); 7];

    % Objective function
    obj_fun = @(vars) COT14(vars(1:end-1), t, vars(end), uField, vField, wField, dt, p, U, Forcing);

    % Simulated annealing options
    options = optimoptions('simulannealbnd', ...
        'Display', 'iter', ...
        'PlotFcns', {@saplotbestf, @saplottemperature, @saplotf, @saplotstopping},...
        'MaxIterations', 1000000, ...
        'MaxFunctionEvaluations', 200000, ...
        'TemperatureFcn', @temperatureboltz, ...
        'InitialTemperature', 5000, ...
        'ReannealInterval', 1000);

    % Perform simulated annealing
    [result, fval] = simulannealbnd(obj_fun, initial_coeffs, lb, ub, options);

    % Extract optimized coefficients and W
    optimized_coeffs = result(1:end-1);
    optimized_W = result(end);

    % Calculate total energy
    totalEnergy = COT14(optimized_coeffs, t, optimized_W, uField, vField, wField, dt, p, U, Forcing);
end
