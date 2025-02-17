% close all;clear,clc;

% Define the file path
file_path = '/Users/boyi/Simulations/hit/output_fields';


% Define the range of file numbers
start_num = 350;
end_num = 1050;

% Initialize cell arrays to store the results
dims_list = cell(end_num - start_num + 1, 1);
valnames_list = cell(end_num - start_num + 1, 1);
values_list = cell(end_num - start_num + 1, 1);
varnames_list = cell(end_num - start_num + 1, 1);
variables_list = cell(end_num - start_num + 1, 1);

% Loop through the range of file numbers
for i = start_num:end_num
    % Construct the filename with leading zeros
    filename = fullfile(file_path, sprintf('field_%05d', i));
    
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
numb=700;
numb1=701;

dt = 1e-3;
t = 0:dt:numb*dt;
W = 4; 
p =1;
tau_p=1e-2;
Forcing = 1;


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
% tic
% 
% 
% %[optimized_coeffs1, optW1, totalEnergy1, fval1] = optimization27(t, W, uField, vField, wField, dt, p, U,Forcing, 10000, 22);  
% 
% 
% [optimized_coeffs2, optW2, totalEnergy2, fval2] = optimization27(t, W, uField, vField, wField, dt, p, U,Forcing, 30000, 22); 
% 
% % 
% % [optimized_coeffs3, optW3, totalEnergy3, fval3] = optimization27(t, W, uField, vField, wField, dt, p, U,Forcing, 10000); 
% % 
% % 
% % [optimized_coeffs4, optW4, totalEnergy4, fval4] = optimization27(t, W, uField, vField, wField, dt, p, U,Forcing,10000); 
% % 
% % [optimized_coeffs5, optW5, totalEnergy5, fval5] = optimization27(t, W, uField, vField, wField, dt, p, U,Forcing, 13000); 
% % 
% % 
% % [optimized_coeffs6, optW6, totalEnergy6, fval6] = optimization27(t, W, uField, vField, wField, dt, p, U,Forcing, 13000);
% 
[optimized_coeffs7, optW7, totalEnergy7, fval7] = optimization27(t(1:300), W, uField(1:300), vField(1:300), wField(1:300), dt, p, U,Forcing, 30000, 12);
% 
% %[optimized_coeffs8, optW8, totalEnergy8, fval8] = optimization27(t(1:300), W, uField(1:300), vField(1:300), wField(1:300), dt, p, U,Forcing, 15000, 12);
% 
% %[optimized_coeffs9, optW9, totalEnergy9, fval9] = optimization27(t(1:300), W, uField(300:600), vField(300:600), wField(300:600), dt, p, U,Forcing, 10000, 12);
% 
% toc


%%


%parpool; % Start parallel pool

N = 14; % Number of parallel searches
results = cell(1, N);
fvals = cell(1, N);

results3h30k = cell(1, N);
fvals3h30k = cell(1, N);

results6 = cell(1, N);
fvals6 = cell(1, N);


G = [0.01;0.05;0.1;0.2;0.4;0.5;0.6;0.7;0.9;1.0;1.2;1.7;2;2.5];

Ghigh = [2,2.5,3,4,5,6,7,8,9,10,11,15,20,25];

parfor i = 1:N
    [results{i}.coeffs, results{i}.W, results{i}.energy, fvals{i}] = ...
        optimization27(t(1:700), W, uField(1:700), vField(1:700), wField(1:700), dt, p, U, Forcing, 30000, 20, G(i),2, 7);
end

parfor i = 1:N
    [results6{i}.coeffs, results6{i}.W, results6{i}.energy, fvals6{i}] = ...
        optimization27(t(1:700), W, uField(1:700), vField(1:700), wField(1:700), dt, p, U, Forcing, 30000, 20, G(i),2,7);
end

parfor i = 1:N
    [results7hhigh{i}.coeffs, results7hhigh{i}.W, results7hhigh{i}.energy, fvals7hhigh{i}] = ...
        optimization27(t(1:700), W, uField(1:700), vField(1:700), wField(1:700), dt, p, U, Forcing, 30000, 20, Ghigh(i),3,15);
end


%%
for i = 3:length(G)
    
    % Plot each point
    plot(G(i), double(fvals1{i})/1.611, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    hold on; % Hold the plot to overlay multiple points
end
set(gca, 'XScale', 'log');
% Adding labels and title
xlabel('G');
ylabel('cot');
title('Plot of cot vs G');
grid on;
hold off; % Release the plot hold



%%
% Given data

% Initialize arrays for cot and G
cot = zeros(1, length(G)); % Preallocate cot array
am = 4;
% Use a for loop to iterate over each value of G
for i = 1:length(G)
    % Call COT14 for each value of G(i)
    then = i;
    [cot(then)] = COT14(results1{am}.coeffs,t(1:300),results1{am}.W, uField(1:300),vField(1:300),wField(1:300),dt,p,U,Forcing, G(then));
    
    % Plot each point
    plot(G(i), cot(i)/1.611, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    hold on; % Hold the plot to overlay multiple points
end

% Adding labels and title
xlabel('G');
ylabel('cot');
title('Plot of cot vs G');
grid on;
hold off; % Release the plot hold







%%
then = 12;
thecoeffs = results1{then}.coeffs;
thew = results1{then}.W;
[cot, fd, gt, ali] = COT14(results1{then}.coeffs,t(1:300),results1{then}.W, uField(1:300),vField(1:300),wField(1:300),dt,p,U,Forcing, G(then));
% [cot, FD, G, ali1] = COT14(optimized_coeffs2,t,optW, uField,vField,wField,dt,p,U,Forcing);
% 
% [cot1, FD, G, ali1] = COT14(optimized_coeffs1k,t,optW1k, uField,vField,wField,dt,p,U,Forcing);
% [cot5, FD, G, ali5] = COT14(optimized_coeffs5k,t,optW5k, uField,vField,wField,dt,p,U,Forcing);
% [cot8, FD, G, ali8] = COT14(optimized_coeffs8k,t,optW8k, uField,vField,wField,dt,p,U,Forcing);

%%
drag = sum((vecnorm(fd).^2).^(3/4));
tail = sum(vecnorm(ali))/length(t(1:300));
%tail1 = sum(vecnorm(ali1))/length(t);

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

    X(:,:) = X14(t, thecoeffs,thew);
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
title('Trajectory, G=0.24');
xlabel('X/L');
ylabel('Y/L');
zlabel('Z/L');
nub=1;
L = 0.17; % Define L for normalization

% Generate colors for the trajectories
colors = lines(5); % Generate 5 distinct colors

% Loop through the first 5 trajectories
for k = 1:1
    % Compute the trajectory for each set
    X(:,:) = X14(t, thecoeffs,thew);
    % Extract X, Y, Z coordinates and normalize by L
    X_k = X(1, :, k) / L;
    Y_k = X(2, :, k) / L;
    Z_k = X(3, :, k) / L;
    
    % Plot the trajectory and assign legend labels
    plot3(X_k, Y_k, Z_k, 'LineWidth', 2, 'Color', colors(k, :), 'DisplayName', ['Trajectory ' num2str(k)]);
    
    % Calculate velocity vectors for the current trajectory
    V = get_vel(uField, vField, wField, 64, X(:, :, k));
    
    % Scale vectors
    V_scale = 1;
    U_ = V(1, :) * V_scale / L; % X-components of vectors
    V_ = V(2, :) * V_scale / L; % Y-components of vectors
    W_ = V(3, :) * V_scale / L; % Z-components of vectors
    
    % Plot vectors at every 10th time step (no DisplayName for vectors)
    step = 10; % Vector plotted every 10th time step
    quiver3(X_k(1:step:end), Y_k(1:step:end), Z_k(1:step:end), ...
            U_(1:step:end), V_(1:step:end), W_(1:step:end), 1, 'Color', colors(k, :), 'LineWidth', 1);
end

% Set specific axes limits (adjusted to account for L)
xlim([0 1/L]);  % Normalized X-axis limits
ylim([0 2/L]);  % Normalized Y-axis limits
zlim([0 1/L]);  % Normalized Z-axis limits

axis equal
% Adjust the view angle for better 3D perception
view(3); % Default 3D view

% Add legend for trajectories only
legend('show', 'Location', 'best');

hold off; % Release the hold on the current figure


% 
% save('resultfeb8first650.mat', 'optimized_coeffs');
% save('Wfeb8first650', 'optW')


%%
xasix = [20000, 5000, 5000, 10000, 10000, 3000, 3000];
yasix = [cot/cot6,cot1/cot6,cot2/cot6,cot3/cot6,cot4/cot6,cot5/cot6,cot6/cot6];
plot(xasix, yasix,'o')
title('Trajectory, 10k trial points');
xlabel('Number of Trialpoints');
ylabel('COT normalized by largest COT');

%%

% Define the filename
filename = 'feb17data.mat';

% Save optimized_coeffs1 to optimized_coeffs9 and optW1 to optW9
save(filename, 'results1', 'results2', 'results3', ...
               'results3h30k', 'results7h30k1', 'results7h30k2');




