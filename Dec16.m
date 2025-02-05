close all;clear,clc;

% Define the file path
file_path = '/home/boyi/Documents/Fields/output_fields/';


% Define the range of file numbers
start_num = 350;
end_num = 3000;

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
numb=300;
numb1=301;

dt = 1e-3;
t = 0:dt:numb*dt;
W = 7; 
p =1;
tau_p=1e-2;
Forcing = 1;

StartLoc = [0.3, 0.5];
coeze = zeros(5,1);

%%
uField = cell(numb, 1);
vField = cell(numb, 1);
wField = cell(numb, 1);

for i = 1:numb1
    data = variables_list{i};
    uField{i} = data(1);
    vField{i} = data(2);
    wField{i} = data(3);
end

% Ufactor = (1/6.975);
% Ufactor = 1;
% 
% 
% [uField, vField, wField] = ChangeU(uField,vField,wField, Ufactor);

U = calculateRMS(uField(1:200),vField(1:200),wField(1:200));
%%
tic
[optimized_coeffs, optW, totalEnergy, fval] = optimization27(t, W, uField, vField, wField, dt, p, U,Forcing);  
toc


%%
[cot, FD, G, ali] = COT14(opt1(1:26),t,opt1(27), uField,vField,wField,dt,p,U,Forcing);

%%

tail = sum(vecnorm(ali))/length(t);

%%

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

    X(:,:) = X14(t, optimized_coeffs,optW);
    plot3(X(1, :), X(2, :), X(3, :), 'LineWidth', 2, 'Color', colors(1, :));



% Set specific axes limits
% xlim([0 1]);  % X-axis limits
% ylim([0 1]);  % Y-axis limits
% zlim([0 1]);  % Z-axis limits

% Adjust the view angle for better 3D perception
view(3); % Default 3D view

hold off; % Release the hold on the current figure


%%
figure;
hold on; % Hold on to plot multiple trajectories in the same figure
grid on; % Enable grid
title('Trajectory, G=0.014');
xlabel('X Component');
ylabel('Y Component');
zlabel('Z Component');
nub=1;
% Generate colors for the trajectories
colors = lines(5); % Generate 5 distinct colors

% Loop through the first 5 trajectories
for k = 1:1
    % Compute the trajectory for each set
    X(:,:,k) = X14(t, optimized_coeffs,optW);
    % Extract X, Y, Z coordinates
    X_k = X(1, :, k);
    Y_k = X(2, :, k);
    Z_k = X(3, :, k);
    
    % Plot the trajectory and assign legend labels
    plot3(X_k, Y_k, Z_k, 'LineWidth', 2, 'Color', colors(k, :), 'DisplayName', ['Trajectory ' num2str(k)]);
    
    % Calculate velocity vectors for the current trajectory
    V = get_vel(uField, vField, wField, 64, X(:, :, k));
    
    % Scale vectors
    V_scale = 1;
    U_ = V(1, :) * V_scale; % X-components of vectors
    V_ = V(2, :) * V_scale; % Y-components of vectors
    W_ = V(3, :) * V_scale; % Z-components of vectors
    
    % Plot vectors at every 10th time step (no DisplayName for vectors)
    step = 10; % Vector plotted every 10th time step
    quiver3(X_k(1:step:end), Y_k(1:step:end), Z_k(1:step:end), ...
            U_(1:step:end), V_(1:step:end), W_(1:step:end), 1, 'Color', colors(k, :), 'LineWidth', 1);
end

% Set specific axes limits
% xlim([0 1]);  % X-axis limits
% ylim([0 1]);  % Y-axis limits
% zlim([0 1]);  % Z-axis limits

% Adjust the view angle for better 3D perception
view(3); % Default 3D view

% Add legend for trajectories only
legend('show', 'Location', 'best');

hold off; % Release the hold on the current figure


save('result10klin.mat', 'optimized_coeffs');
save('W10klin', 'optW')



