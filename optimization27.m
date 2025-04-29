function [optimized_coeffs, optimized_W, totalEnergy, fval] = optimization27(t, W, uField, vField, wField, dt, p, U, Forcing, trialpoints, mode, tau_p, wlow, whigh, lim, ip, prev_coeffs)

    % mode = mode;
    

% Default zero init
initial_coeffs = zeros(3*mode+3, 1);
initial_coeffs(end) = W;

if exist('prev_coeffs', 'var') && ~isempty(prev_coeffs)
    % Total length minus 2 control inputs (no W)
    prev_mode = (length(prev_coeffs) - 2) / 3;

    % Check that prev_mode is valid
    if mod(prev_mode, 1) == 0 && prev_mode <= mode
        % Copy X coefficients
        initial_coeffs(1:prev_mode) = prev_coeffs(1:prev_mode);

        % Copy Y coefficients
        initial_coeffs(mode+1:mode+prev_mode) = prev_coeffs(prev_mode+1:2*prev_mode);

        % Copy Z coefficients
        initial_coeffs(2*mode+1:2*mode+prev_mode) = prev_coeffs(2*prev_mode+1:3*prev_mode);

        % Copy control terms
        initial_coeffs(end-2:end-1) = prev_coeffs(end-1:end);  % last two elements of prev_coeffs
    else
        warning('[Init] prev_coeffs length mismatch or invalid mode!');
    end
end




    
    % Define bounds
    lb = [-lim * ones(mode, 1); -lim * ones(mode, 1); -lim * ones(mode, 1); -ip*ones(2,1) ;wlow];
    ub = [lim * ones(mode, 1); lim * ones(mode, 1); lim * ones(mode, 1); ip*ones(2,1) ;whigh];



    % Define optimization options using the interior-point algorithm
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'none');

    % Define the optimization problem for fmincon
    problem = createOptimProblem('fmincon', 'objective', ...
        @(vars) COT14(vars(1:end-1), t, vars(end), uField, vField, wField, dt, p, U, Forcing, tau_p), ...
        'x0', initial_coeffs, 'lb', lb, 'ub', ub, 'options', options);

    % Create a GlobalSearch object to perform the optimization
      gs = GlobalSearch('NumTrialPoints', trialpoints);


    % === Run GlobalSearch optimization ===
[result_gs, fval_gs] = run(gs, problem);
coeffs_gs = result_gs(1:end-1);
W_gs = result_gs(end);
energy_gs = COT14(coeffs_gs, t, W_gs, uField, vField, wField, dt, p, U, Forcing, tau_p);

% === Evaluate warm-start point directly (prev_coeffs padded into initial_coeffs) ===
energy_warm = Inf;
if exist('prev_coeffs', 'var') && ~isempty(prev_coeffs)
    coeffs_warm = initial_coeffs(1:end-1);
    W_warm = initial_coeffs(end);
    energy_warm = COT14(coeffs_warm, t, W_warm, uField, vField, wField, dt, p, U, Forcing, tau_p);
end

% === Compare and select the best result ===
[~, winner] = min([energy_warm, energy_gs]);

if winner == 1
    optimized_coeffs = coeffs_warm;
    optimized_W = W_warm;
    totalEnergy = energy_warm;
    fval = energy_warm;  % Not really fval, but acceptable for tracking
    % fprintf('[Use] Warm-start (padded) solution: %.4f\n', energy_warm);
else
    optimized_coeffs = coeffs_gs;
    optimized_W = W_gs;
    totalEnergy = energy_gs;
    fval = fval_gs;
    % fprintf('[Use] GlobalSearch result: %.4f\n', energy_gs);
end


    % Display the result
    % disp('Optimized total energy:');
    % disp(totalEnergy);
end

