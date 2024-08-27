close all; clc; clear;
%the flow flow is a 1*1*1 box
%set tf to be 1 and 1000 time step.
t = linspace(0,1, 1000);
%velocity in y is 1
W = 0.7;

StartLoc = [0.5, 0.5];


% Define the coefficients and parameters
% Number of sets
num_sets = 1000;

% Define the magnitude coefficient
coe = 0.1;  % Adjust this to change the range of randomness

% Initialize B with zeros including first and third index for all sets
B = zeros(4, num_sets);

% Generate random numbers for A with range [-coe, coe]
A = (2 * coe) * rand(4, num_sets) - coe;

% Generate random numbers for C with range [-coe, coe]
C = (2 * coe) * rand(4, num_sets) - coe;

% Generate random numbers for the second and fourth indices of B
B(2, :) = (2 * coe) * rand(1, num_sets) - coe;
B(4, :) = (2 * coe) * rand(1, num_sets) - coe;
tf = 1;

% Initialize a 4D array to hold the results: 3 dimensions for x, variable length of t, 1000 sets
X = zeros(3, length(t), 1000);

% Loop over the 1000 sets of coefficients
for k = 1:1000
    % Extract the k-th set of coefficients
    a_k = A(:, k);
    b_k = B(:, k);
    c_k = C(:, k);

    % Initialize temporary storage for x(t) for the current set of coefficients
    x_temp = zeros(3, length(t));
    
    % Compute x(t) for each time point
    for i = 1:length(t)
        x_temp(:, i) = calculateX(t(i), a_k, b_k, c_k, tf, W, StartLoc(1), StartLoc(2));
    end

    % Store the result in the 4D array
    X(:, :, k) = x_temp;
end

% Define the figure for plotting
figure;
hold on; % Hold on to plot multiple trajectories in the same figure
grid on; % Enable grid
title('3D Trajectory of x(t) for First 10 Sets');
xlabel('X_1 Component');
ylabel('X_2 Component');
zlabel('X_3 Component');

% Loop through the first 10 trajectories
colors = lines(10); % Generate 10 distinct colors
for k = 1:10
    plot3(X(1, :, k), X(2, :, k), X(3, :, k), 'LineWidth', 2, 'Color', colors(k, :));
end

% Set specific axes limits
xlim([0 1]);  % X-axis limits
ylim([0 1]);  % Y-axis limits
zlim([0 1]);  % Z-axis limits

% Adjust the view angle for better 3D perception
view(3); % Default 3D view

hold off; % Release the hold on the current figure





