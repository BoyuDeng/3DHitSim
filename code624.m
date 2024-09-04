close all; clear; clc;
file_path = '/Users/boyi/Research/Research/field/';


% Define the range of file numbers
start_num = 350;
end_num = 1000;

% Initialize cell arrays to store the results
dims_list = cell(end_num - start_num + 1, 1);
valnames_list = cell(end_num - start_num + 1, 1);
values_list = cell(end_num - start_num + 1, 1);
varnames_list = cell(end_num - start_num + 1, 1);
variables_list = cell(end_num - start_num + 1, 1);

% Loop through the range of file numbers
for i = start_num:end_num
    % Construct the filename with leading zeros
    filename = sprintf('field_%05d', i);
    
    % Read the file
    [dims, valnames, values, varnames, variables] = read_field_file(filename);
    
    % Store the results
    dims_list{i - start_num + 1} = dims;
    valnames_list{i - start_num + 1} = valnames;
    values_list{i - start_num + 1} = values;
    varnames_list{i - start_num + 1} = varnames;
    variables_list{i - start_num + 1} = variables;
end

%%
uField = cell(651, 1);
vField = cell(651, 1);
wField = cell(651, 1);

for i = 1:651
    data = variables_list{i};
    uField{i} = data(1);
    vField{i} = data(2);
    wField{i} = data(3);
end

U = calculateRMS(uField,vField,wField);

%%


dt = 1e-3;
%the flow flow is a 1*1*1 box
%set tf to be 1 and 1000 time step.
t = 0:dt:650*dt;
%velocity in y is 2
W = 1; 
p =2;
Forcing = 0;

StartLoc = [0.5, 0.5];

coe = [0,0,0];

[straightCOT, straightFdrag] = COT_function(coe,t,W,StartLoc,uField,vField,wField,dt, p,U, Forcing);

%%

% Define the coefficients and parameters
% Number of sets
num_sets = 1000;

% Define the magnitude coefficient
coe = 0.1;  % Adjust this to change the range of randomness


% Generate random numbers for A with range [-coe, coe]
A = (2 * coe) * rand(12, num_sets) - coe;


% Define the figure for plotting
figure;
hold on; % Hold on to plot multiple trajectories in the same figure
grid on; % Enable grid
title('3D Trajectory of x(t) for First 5 Sets');
xlabel('X Component');
ylabel('Y Component');
zlabel('Z Component');

% Loop through the first 10 trajectories
colors = lines(10); % Generate 10 distinct colors
for k = 1:5
    X(:,:,k) = calculateX(t, A(:,k),W, StartLoc);
    plot3(X(1, :, k), X(2, :, k), X(3, :, k), 'LineWidth', 2, 'Color', colors(k, :));
end



% Set specific axes limits
xlim([0 1]);  % X-axis limits
ylim([0 1]);  % Y-axis limits
zlim([0 1]);  % Z-axis limits

% Adjust the view angle for better 3D perception
view(3); % Default 3D view

hold off; % Release the hold on the current figure



%%
E = zeros(1,1000);

for i = 1:1000
    E(i) = COT_function(A(:,i),t,W,StartLoc,uField,vField,wField,dt,p,U,Forcing);
end

histogram(E/straightCOT)
title('Histogram of 1000 trajectories');
xlabel('Normalized COT');
ylabel('Frequency');
%%
random_numbers = rand(1, 100);

% Scale and shift to the range 0.3 to 0.7
scaled_random_numbers = 0.3 + (0.7 - 0.3) * random_numbers;
random_numbers1 = rand(1, 100);

% Scale and shift to the range 0.3 to 0.7
scaled_random_numbers1 = 0.3 + (0.7 - 0.3) * random_numbers1;

StartLoc = [scaled_random_numbers; scaled_random_numbers1];
coeze = zeros(3,1);
for i = 1:100
    Es(i) = COT_function(coeze,t,W,StartLoc(:,i),uField,vField,wField,dt,p,U,Forcing);
end
histogram(Es/straightCOT)
title('Histogram of 100 strating point');
xlabel('Normalized COT');
ylabel('Frequency');
%%
[min_value, min_index] = min(E);
% Define bounds for the coefficients
% lb = -0.1 * ones(12, 1);  % Lower bounds
% ub = 0.1*ones(12, 1);   % Upper bounds
% initial_coeffs = A(:,min_index);
% % Optimization options
% options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter');
% 
% % Perform optimization using fmincon
% result = fmincon(@(coeffs) COT_function(coeffs, t, W, StartLoc, uField, vField, wField, dt, p), ...
%                  initial_coeffs, [], [], [], [], lb, ub, [], options);
% 
% % Optimized coefficients
% optimized_coeffs = result;
% 
% 
% totalEnergy = COT_function(optimized_coeffs, t, W, StartLoc, uField, vField, wField, dt, p);
% disp('Initial total energy:');
% disp(totalEnergy);

%%
% Define the problem for fmincon
initial_coeffs = A(:,min_index);
lb = -0.1 * ones(12, 1);  % Lower bounds
ub = 0.1*ones(12, 1);   % Upper bounds
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter');
problem = createOptimProblem('fmincon', 'objective', ...
    @(coeffs) COT_function(coeffs, t, W, StartLoc, uField, vField, wField, dt, p, U,Forcing), ...
    'x0', initial_coeffs, 'lb', lb, 'ub', ub, 'options', options);

% Create a GlobalSearch object
gs = GlobalSearch;

% Run the optimization
[result, fval] = run(gs, problem);

% Optimized coefficients
optimized_coeffs = result;

% Calculate the total energy with the optimized coefficients
totalEnergy = COT_function(optimized_coeffs, t, W, StartLoc, uField, vField, wField, dt, p, U,Forcing);
disp('Optimized total energy:');
disp(totalEnergy);

%%
facotor1 = 3;
% Preallocate cells to store the multiplied data
    uFieldMultiplied = cell(size(uField));
    vFieldMultiplied = cell(size(vField));
    wFieldMultiplied = cell(size(wField));

    % Convert cell arrays to matrices using a loop and multiply by factor
    for t = 1:length(t)
        u = cell2mat(uField{t}) .* facotor1;
        v = cell2mat(vField{t}) .* facotor1;
        w = cell2mat(wField{t}) .* facotor1;

        if isnumeric(u) && isnumeric(v) && isnumeric(w)
            uFieldMultiplied{t} = mat2cell(u, size(u,1), size(u,2), size(u,3));
            vFieldMultiplied{t} = mat2cell(v, size(v,1), size(v,2), size(v,3));
            wFieldMultiplied{t} = mat2cell(w, size(w,1), size(w,2), size(w,3));
        else
            error('Each cell must contain numeric data.');
        end
    end

    %%

    % Normalize the results by straightCOT
normalized_E = E / straightCOT;
normalized_E_distinct = totalEnergy / straightCOT;

% Plot the histogram for the 1000 trajectories
figure;
histogram(normalized_E);
hold on;

% Mark the distinct point on the histogram as a vertical line
xline(normalized_E_distinct, 'r', 'LineWidth', 2);

% Add title and labels
title('Histogram of 1000 trajectories with distinct point');
xlabel('Normalized COT');
ylabel('Frequency');

% Add legend for the distinct point
legend('1000 Trajectories', 'Distinct Point', 'Location', 'Best');
hold off;




