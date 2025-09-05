% Define the objective function (Rosenbrock function)
objectiveFunction = @(x) (1 - x(1))^2 + 100 * (x(2) - x(1)^2)^2;

% Define the bounds for the variables
lb = [-5, -5]; % Lower bounds
ub = [5, 5];   % Upper bounds

% Define initial guess
x0 = [0, 0]; % Starting point for the algorithm

% Set options for simulated annealing
options = optimoptions('simulannealbnd', ...
    'Display', 'iter', ...
    'PlotFcns', {@saplotbestf, @saplottemperature, @saplotf, @saplotstopping});

% Run the simulated annealing algorithm
[xOpt, fVal] = simulannealbnd(objectiveFunction, x0, lb, ub, options);

% Display results
disp('Optimal solution found:');
disp(xOpt);
disp('Function value at optimal solution:');
disp(fVal);

