close all;clear,clc;

load("fielddata624.mat");

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

Ufactor = (1/6.975)*0.0098;

[uField, vField, wField] = ChangeU(uField,vField,wField, Ufactor);


U = calculateRMS(uField,vField,wField);



%%

dt = 1e-3;
%the flow flow is a 1*1*1 box
%set tf to be 1 and 1000 time step.
t = 0:dt:650*dt;
%velocity in y is 2
W = 1; 
p =1;


StartLoc = [0.5, 0.5];
coeze = zeros(5,1);


%%


% Define the coefficients and parameters
% Number of sets
num_sets = 1000;

% Define the magnitude coefficient
coe = 0.1;  % Adjust this to change the range of randomness


% Generate random numbers for A with range [-coe, coe]
A = (2 * coe) * rand(14, num_sets) - coe;


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
    X(:,:,k) = X14(t, A(:,k),W);
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
FD = zeros(1,1000);
Forcing = 1;
for i = 1:1000
    E(i) = COT14(A(:,i),t,W,uField,vField,wField,dt,p,U,Forcing);
end
[straightCOT, G] = COT14(coeze,t,W,uField,vField,wField,dt, p,U, Forcing);
histogram(E/straightCOT)
title('Histogram of 1000 random trajectories');
xlabel('Normalized COT');
ylabel('Frequency');


%%

[optimized_coeffs, optimized_E] = optimization14(t, W, uField, vField, wField, dt, p, U,Forcing);


%%

    % Normalize the results by straightCOT
normalized_E = E / straightCOT;
normalized_E_distinct = optimized_E / straightCOT;

% Plot the histogram for the 1000 trajectories
figure;
histogram(normalized_E);
hold on;

% Mark the distinct point on the histogram as a vertical line
xline(normalized_E_distinct, 'r', 'LineWidth', 2);

% Add title and labels
title('Histogram of 1000 trajectories with optimized point');
xlabel('Normalized COT');
ylabel('Frequency');

% Add legend for the distinct point
legend('1000 Trajectories', 'Optimal Trajecotry', 'Location', 'Best');
hold off;


%%
% Example variables
% P is the 3xT matrix for particle positions
% V is the 3xT matrix for vectors corresponding to the particle

t = 0:dt:650*dt;
P = X14(t,optimized_coeffs,W);
V = get_vel(uField,vField,wField,64,P);


% Assume P and V are already defined
t = size(P, 2); % number of time steps

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
        U_(1:step:end), V_(1:step:end), W_(1:step:end), 0.005, 'r'); % scale factor 0.5 for the vectors

% Add labels and title
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Plot of Particle Motion and Vectors');
grid on;

% Adjust view angle for better visualization
view(3); % standard 3D view

hold off;





