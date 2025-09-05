close all;clear,clc;

load("fielddata624.mat");
dt = 1e-3;
t = 0:dt:650*dt;
W = 0.6; 
p =1;
tau_p=1e-2;
Forcing = 1;

StartLoc = [0.5, 0.5];
coeze = zeros(5,1);
% Initialize 2D cell arrays with dimensions [651, n]
n = 10; % assuming 'n' is defined somewhere above as the length of Ufactor
uField = cell(651, n);
vField = cell(651, n);
wField = cell(651, n);
orguField = cell(651,1);
orgwField = orguField;
orgvField = orgwField;
G = zeros(n,1);
U = G;
% Populate initial values for the first column
for i = 1:651
    data = variables_list{i};
    orguField{i} = data(1);
    orgvField{i} = data(2);
    orgwField{i} = data(3);
end

% Define Ufactor with n values
Ufactor = linspace(0.005, 0.1, n);

% Call ChangeU function for each Ufactor to fill each column of uField, vField, wField
for j = 1:n
    % Extract the current column (original values) for modification
    current_uField = orguField;
    current_vField = orgvField;
    current_wField = orgwField;
    
    % Apply ChangeU to get the modified fields with current Ufactor(j)
    [updated_uField, updated_vField, updated_wField] = ChangeU(current_uField, current_vField, current_wField, Ufactor(j));
    
    % Store the results in the j-th column of the output fields
    uField(:, j) = updated_uField;
    vField(:, j) = updated_vField;
    wField(:, j) = updated_wField;
    U(j) = calculateRMS(updated_uField,updated_vField,updated_wField);
    G(j) = calculateG(tau_p, U(j));
end


%%
i = 1;
[optimized_coeffs, optimized_W,optimized_E] = optimization27(t, W, uField(:,i), vField(:,i), wField(:,i), dt, p, U(i),Forcing);



%%
for i = 1:10


[optimized_coeffs, optimized_W,optimized_E] = optimization27(t, W, uField(:,i), vField(:,i), wField(:,i), dt, p, U(i),Forcing);
optimized_Es(i) = optimized_E;
optcoeffs(:,i) = optimized_coeffs;
optW(i) = optimized_W;

[straightCOT] = COT14(coeze,t,optimized_W,uField(:,i),vField(:,i),wField(:,i),dt, p,U, Forcing);
straightCOTs(i) = straightCOT;

end

%%

for i =1:10
    G(i) = calculateG(tau_p, U(i));
    [straightCOT] = COT14(coeze,t,optW(i),uField(:,i),vField(:,i),wField(:,i),dt, p,U(i), Forcing);
    straightCOTs(i) = straightCOT;
    result(i) = optimized_Es(i)/straightCOTs(i);
end

%%
figure;
loglog(G, result, '-o'); % '-o' specifies a line with circle markers
xlabel('G');
ylabel('OptCOT/StraightCOT');
title('Plot of ten Gs');
grid on;



%%
save('27var.mat', 'G', 'result', 'optcoeffs','optW','U')


%%

% % Example variables
% P is the 3xT matrix for particle positions
% V is the 3xT matrix for vectors corresponding to the particle
i = 2;
t = 0:dt:650*dt;
P = X14(t,optcoeffs(:,i),W);
V = get_vel(uField(:,i),vField(:,i),wField(:,i),64,P);


% Assume P and V are already defined
%t = size(P, 2); % number of time steps

% Extract X, Y, Z coordinates for the particle locations
X = P(1, :);  % X-coordinates of particle
Y = P(2, :);  % Y-coordinates of particle
Z = P(3, :);  % Z-coordinates of particle


V_scale = 1;
% Extract vector components at each time step
U_ = V(1, :)*V_scale;  % X-components of the vectors
V_ = V(2, :)*V_scale; % Y-components of the vectors
W_ = V(3, :)*V_scale;  % Z-components of the vectors

% Create the figure
figure;

% Plot the particle location in 3D
plot3(X, Y, Z, 'LineWidth', 1.5);
hold on;

% Plot vectors at every 10th time step
step = 10; % Vector plotted every 10th time step
quiver3(X(1:step:end), Y(1:step:end), Z(1:step:end), ...
        U_(1:step:end), V_(1:step:end), W_(1:step:end), 0.5, 'r'); % scale factor 0.5 for the vectors

% Add labels and title
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Plot of Particle Motion and Vectors');
grid on;

% Adjust view angle for better visualization
view(3); % standard 3D view

hold off;

%%

% Create a 3D grid for the original 64x64x64 matrix
[x, y, z] = ndgrid(1:64, 1:64, 1:64);

% Downsample the data for visualization
step = 2; % Adjust this step size based on how much you want to downsample
x = x(1:step:64, 1:step:64, 1:step:64) / 64; % Normalize to [0, 1]
y = y(1:step:64, 1:step:64, 1:step:64) / 64; % Normalize to [0, 1]
z = z(1:step:64, 1:step:64, 1:step:64) / 64; % Normalize to [0, 1]

% Access specific u, v, w fields from cell arrays
nu = 1; % Assuming you are interested in the first set of velocity fields
u = uField{nu, 1}{1}; % Extract the matrix from the cell array
v = vField{nu, 1}{1};
w = wField{nu, 1}{1};

% Downsample the velocity matrices
u = u(1:step:64, 1:step:64, 1:step:64);
v = v(1:step:64, 1:step:64, 1:step:64);
w = w(1:step:64, 1:step:64, 1:step:64);

% Select a 2D slice in the XZ plane (for example, taking y=1 slice)
slice_y = 1;
x_slice = squeeze(x(:, slice_y, :));
z_slice = squeeze(z(:, slice_y, :));
u_slice = squeeze(u(:, slice_y, :));
w_slice = squeeze(w(:, slice_y, :));

% Plot the 2D vector field using quiver
figure;
quiver(x_slice, z_slice, u_slice, w_slice);
xlabel('X');
ylabel('Z');
title('Starting Location, a XZ clip(travel in y direction)');
grid on;
axis([0 1 0 1]); % Set axis limits to [0, 1] for both X and Z

% Define and add specific solid points
hold on;
plot(0.5, 0.5, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'result 1 Center'); % Center point
plot(0.3, 0.3, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'result 2 (0.3, 0.3)');
plot(0.3, 0.7, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'result 3 (0.3, 0.7)');
plot(0.7, 0.7, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'DisplayName', 'result 4 (0.7, 0.7)');
plot(0.7, 0.3, 'co', 'MarkerSize', 8, 'MarkerFaceColor', 'c', 'DisplayName', 'result 5 (0.7, 0.3)');
hold off;

legend('show'); % Display legend for the points


%%
% Transpose optW
optW_transpose = optW';

% Perform element-wise division
result = optW_transpose ./ U;

% Plot G versus the result using dots as markers
loglog(G, result, 'o');

% Add labels and title (optional)
xlabel('G');
ylabel('optW\/U');
title('Plot of G vs optW/U');

% Show grid (optional)
grid on;

%%
figure;
hold on; % Hold on to plot multiple trajectories in the same figure
grid on; % Enable grid
title('Trajectory for G=2.8,0.9,0.53,0.38,0.29');
xlabel('X Component');
ylabel('Y Component');
zlabel('Z Component');

% Generate colors for the trajectories
colors = lines(5); % Generate 5 distinct colors

% Loop through the first 5 trajectories
for k = 1:5
    % Compute the trajectory for each set
    X1(:, :, k) = X14(t, optcoeffs(:, k), optW(k));
    
    % Extract X, Y, Z coordinates
    X_k = X1(1, :, k);
    Y_k = X1(2, :, k);
    Z_k = X1(3, :, k);
    
    % Plot the trajectory and assign legend labels
    plot3(X_k, Y_k, Z_k, 'LineWidth', 2, 'Color', colors(k, :), 'DisplayName', ['Trajectory ' num2str(k)]);
    
    % Calculate velocity vectors for the current trajectory
    V = get_vel(uField(:, 5), vField(:, 5), wField(:, 5), 64, X1(:, :, k));
    
    % Scale vectors
    V_scale = 1;
    U_ = V(1, :) * V_scale; % X-components of vectors
    V_ = V(2, :) * V_scale; % Y-components of vectors
    W_ = V(3, :) * V_scale; % Z-components of vectors
    
    % Plot vectors at every 10th time step (no DisplayName for vectors)
    step = 10; % Vector plotted every 10th time step
    quiver3(X_k(1:step:end), Y_k(1:step:end), Z_k(1:step:end), ...
            U_(1:step:end), V_(1:step:end), W_(1:step:end), 0.5, 'Color', colors(k, :), 'LineWidth', 1);
end

% Set specific axes limits
xlim([0 1]);  % X-axis limits
ylim([0 1]);  % Y-axis limits
zlim([0 1]);  % Z-axis limits

% Adjust the view angle for better 3D perception
view(3); % Default 3D view

% Add legend for trajectories only
legend('show', 'Location', 'best');

hold off; % Release the hold on the current figure




