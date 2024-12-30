close all;clear,clc;

load("fielddata624.mat");

dt = 1e-3;
t = 0:dt:650*dt;
W = 1; 
p =1;
tau_p=1e-2;
Forcing = 1;

StartLoc = [0.3, 0.5];
coeze = zeros(5,1);

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




% Assuming uField and vField are populated with 64x64x64 matrices in each cell
% Initialize parameters
numTimeSteps = 651;   % Total time steps
downsampleFactor = 4; % Downsample for vector density

% Create a figure for the animation
figure;
hold on;

% Pre-allocate the animation frames
frames(numTimeSteps) = struct('cdata', [], 'colormap', []);

% Loop through time steps to create animation frames
for t = 1:numTimeSteps
    % Extract 2D slices at z=1 for u and v, and convert from cell to numeric array
    dataU = uField{t};
    dataV = vField{t};
    
    % Ensure we extract numeric data
    if iscell(dataU)
        dataU = dataU{1}; % Unpack the inner cell
    end
    if iscell(dataV)
        dataV = dataV{1}; % Unpack the inner cell
    end
    
    % Extract the 2D slice for z=1
    u2D = squeeze(dataU(:,:,1));
    v2D = squeeze(dataV(:,:,1));
    
    % Define the grid for the vectors (optional: downsample for clarity)
    [X, Y] = meshgrid(1:size(u2D, 2), 1:size(u2D, 1));
    X = X(1:downsampleFactor:end, 1:downsampleFactor:end);
    Y = Y(1:downsampleFactor:end, 1:downsampleFactor:end);
    u2D = u2D(1:downsampleFactor:end, 1:downsampleFactor:end);
    v2D = v2D(1:downsampleFactor:end, 1:downsampleFactor:end);
    
    % Plot the vector field using quiver
    quiver(X, Y, u2D, v2D, 'k'); % 'k' sets vector color to black
    axis tight;
    axis off;
    title(['Time Step: ', num2str(t)]);
    
    % Capture the current frame
    frames(t) = getframe(gcf);
    
    % Clear the figure for the next frame
    clf;
end

% Play the animation
movie(gcf, frames, 1, 10); % Play animation at 10 frames per second (fps)

% Optionally save as a video
video = VideoWriter('velocity_field_evolution.avi');
video.FrameRate = 10; % Set frame rate
open(video);
writeVideo(video, frames);
close(video);

disp('Animation saved as velocity_field_evolution.avi');

%%
%%
dt=0.1;
t = 0:dt:650*dt;
coe = ones(1, 26);
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

    X(:,:) = X14(t, coe,0.11);
    plot3(X(1, :), X(2, :), X(3, :), 'LineWidth', 2, 'Color', colors(1, :));



% Set specific axes limits


% Adjust the view angle for better 3D perception
view(3); % Default 3D view

hold off; % Release the hold on the current figure



