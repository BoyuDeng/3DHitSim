close all;clear,clc;

load("fielddata624.mat");
dt = 1e-3;
t = 0:dt:650*dt;
W = 1; 
p =1;
tau_p=1e-1;
Forcing = 1;

StartLoc = [0.5, 0.5];
coeze = zeros(5,1);
% Initialize 2D cell arrays with dimensions [651, n]
n = 10; % assuming 'n' is defined somewhere above as the length of Ufactor
uField = cell(651, n);
vField = cell(651, n);
wField = cell(651, n);
G = zeros(n,1);
U = G;
% Populate initial values for the first column
for i = 1:651
    data = variables_list{i};
    uField{i, 1} = data(1);
    vField{i, 1} = data(2);
    wField{i, 1} = data(3);
end

% Define Ufactor with n values
Ufactor = linspace(0.05, 10, n);

% Call ChangeU function for each Ufactor to fill each column of uField, vField, wField
for j = 1:n
    % Extract the current column (original values) for modification
    current_uField = uField(:, 1);
    current_vField = vField(:, 1);
    current_wField = wField(:, 1);
    
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
for i = 1:10


[optimized_coeffs, optimized_W,optimized_E] = optimization15(t, W, uField(:,i), vField(:,i), wField(:,i), dt, p, U(i),Forcing);
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

% Example variables
% P is the 3xT matrix for particle positions
% V is the 3xT matrix for vectors corresponding to the particle

t = 0:dt:650*dt;
P = X14(t,optcoeffs(:,3),W);
V = get_vel(uField(:,3),vField(:,3),wField(:,3),64,P);


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






