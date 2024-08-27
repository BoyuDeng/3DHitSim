
% file_path = '/home/boyi/Desktop/CFD/field/';
close all; clear; clc;
load('newfielddata.mat')
% Define the range of file numbers
% start_num = 350;
% end_num = 1000;
% 
% % Initialize cell arrays to store the results
% dims_list = cell(end_num - start_num + 1, 1);
% valnames_list = cell(end_num - start_num + 1, 1);
% values_list = cell(end_num - start_num + 1, 1);
% varnames_list = cell(end_num - start_num + 1, 1);
% variables_list = cell(end_num - start_num + 1, 1);
% 
% % Loop through the range of file numbers
% for i = start_num:end_num
%     % Construct the filename with leading zeros
%     filename = sprintf('field_%05d', i);
% 
%     % Read the file
%     [dims, valnames, values, varnames, variables] = read_field_file(filename);
% 
%     % Store the results
%     dims_list{i - start_num + 1} = dims;
%     valnames_list{i - start_num + 1} = valnames;
%     values_list{i - start_num + 1} = values;
%     varnames_list{i - start_num + 1} = varnames;
%     variables_list{i - start_num + 1} = variables;
% end
% 
% % Display a message indicating that the files have been read
% disp('All files have been read successfully.');

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


%%


% Create a grid that matches the dimensions of your velocity matrices
[x, y, z] = ndgrid(1:64, 1:64, 1:64);

% Downsample the data for visualization if needed (optional)
% This is useful if the dataset is too large to visualize clearly
step = 1; % Example step size for downsampling
x = x(1:step:10, 1:step:10, 1:step:10);
y = y(1:step:10, 1:step:10, 1:step:10);
z = z(1:step:10, 1:step:10, 1:step:10);
nu=515;
% Assume uField, vField, and wField are cell arrays of matrices
u = uField{nu};
u = u{1};
u = u(1:step:10, 1:step:10, 1:step:10);
v = vField{nu};
v = v{1};
v = v(1:step:10, 1:step:10, 1:step:10);
w = wField{nu};
w = w{1};
w = w(1:step:10, 1:step:10, 1:step:10);


% Create the 3D vector field plot
figure;
quiver3(x, y, z, u, v, w);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Vector Field');
grid on;
axis tight;
view(3); % Set the view to 3D

%%


nu1 = 1;
nu600 = 600;
step = 1; % Example step size for downsampling

% Extract and downsample the data for nu = 1
u1 = uField{nu1};
u1 = u1{1};
u1 = u1(1:step:10, 1:step:10, 1:step:10);
v1 = vField{nu1};
v1 = v1{1};
v1 = v1(1:step:10, 1:step:10, 1:step:10);
w1 = wField{nu1};
w1 = w1{1};
w1 = w1(1:step:10, 1:step:10, 1:step:10);

% Extract and downsample the data for nu = 600
u600 = uField{nu600};
u600 = u600{1};
u600 = u600(1:step:10, 1:step:10, 1:step:10);
v600 = vField{nu600};
v600 = v600{1};
v600 = v600(1:step:10, 1:step:10, 1:step:10);
w600 = wField{nu600};
w600 = w600{1};
w600 = w600(1:step:10, 1:step:10, 1:step:10);

% Calculate the difference field
u_diff = u600 - u1;
v_diff = v600 - v1;
w_diff = w600 - w1;

% Create the grid that matches the dimensions of your velocity matrices
[x, y, z] = ndgrid(1:10, 1:10, 1:10);

% Create the 3D vector field plot of the difference
figure;
quiver3(x, y, z, u_diff, v_diff, w_diff);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Difference of 3D Vector Fields (nu = 600 - nu = 1)');
grid on;
axis


%%

% t_index = 100; % for example, the 100th time step
% i = 32;
% j = 32;
% k = 32;
% 
% % Get the vector at the specified time step and location
% vector = get_vector_at_location(uField, vField, wField, t_index, i, j, k);


dt = 5e-4;
%the flow flow is a 1*1*1 box
%set tf to be 1 and 1000 time step.
t = 0:dt:650*dt;
%velocity in y is 1
W = 2;

StartLoc = [0.5, 0.5];

xstraight = calculateXold(t, 0,0,0,t(end),W, StartLoc(1),StartLoc(2),0);

%xstraight = X(:,:,2);

straightloc = get_particle_cell_Location(xstraight,64);

vectors = get_vectors_at_locations(uField, vField, wField, straightloc);

location_differences = diff(xstraight, 1, 2);
Vo = location_differences/dt;
coe = [0,0,0];
Vo = calculateV(t,coe,W);
Fd = zeros(3,650);
for i = 1:650
    Fd(:,i) = (vectors(:,i)-Vo(:,i))*norm(vectors(:,i)-Vo(:,i));
    Energy(i) = dot(Fd(:,i),location_differences(:,i));
end

% vectors1 = vectors(:,1:end-1);
% Fdrag = (vectors1 - Vo) .* vecnorm(vectors1 - Vo);
% Totalenergy1 = dot(Fdrag, location_differences);
% Totalenergy2 = sum(Totalenergy1);
        

Energy1 = sum(Energy);


% Calculate differences between successive columns


% Calculate velocity (change in location divided by time step)



Realdrag = Energy1/(0.650);
%%

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
        x_temp(:, i) = calculateXold(t(i), a_k, b_k, c_k, t(end), W, StartLoc(1), StartLoc(2), 0.1);
    end

    % Store the result in the 4D array
    X(:, :, k) = x_temp;
end


%%
Particles_locations = zeros(3, 651, length(X));
Velocity_fields = Particles_locations;
Fdrag = zeros(3,650,length(X));
Totalenergy1 = zeros(650,1000);
Totalenergy = zeros(1000,1);
Realdrags = Totalenergy;
Gap = zeros(3,650,length(X));
Vs = Gap;
for i = 1:length(X)
Particles_locations(:,:,i) = get_particle_cell_Location(X(:,:,i),64);
Velocity_fields(:,:,i) = get_vectors_at_locations(uField, vField, wField,Particles_locations(:,:,i));
Gap(:,:,i) = diff(X(:,:,i), 1, 2);
Vs(:,:,i) = Gap(:,:,i)/dt;
    for j = 1:650
        Fdrag(:,j,i)=(Velocity_fields(:,j,i)-Vs(:,j,i))*norm(Velocity_fields(:,j,i)-Vs(:,j,i));
        Totalenergy1(j,i) = dot(Fdrag(:,j,i), Gap(:,j,i));
    end

    Totalenergy(i) = sum(Totalenergy1(:,i));
    Realdrags(i)=Totalenergy(i)/(0.650);


end

NormDrag = Realdrags/Realdrag;

histogram(NormDrag);
title('Histogram of 1000 trajectories');
xlabel('Normalized COT');
ylabel('Frequency');
grid on;  % Add a grid for better readability

% Additional customization (optional)
% histogram(data, 'BinWidth', 0.5);  % Adjust bin width
% histogram(data, 'FaceColor', 'red');  % Change bar color
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
for k = 1:5
    plot3(X(1, :, k), X(2, :, k), X(3, :, k), 'LineWidth', 2, 'Color', colors(k, :));
end



% Set specific axes limits
xlim([0 1]);  % X-axis limits
ylim([0 1]);  % Y-axis limits
zlim([0 1]);  % Z-axis limits

% Adjust the view angle for better 3D perception
view(3); % Default 3D view

hold off; % Release the hold on the current figure


q1 = Gap(:,:,10)
q2 = X(:,:,10)

for q = 1:1000
    for w = 1:650
    dis(q,w) = norm(Gap(:,w,q));
    end
end

ans1 = mean(sum(dis,2))



Particles_locations(:,:,3)


for k = 1:1000
for c = 1:650
Fdrag1(:,c,k) = Fdrag(:,c,k)./Fd(:,c);
end
end

histogram(Fdrag1)
max(max(max(Fdrag1)))

for i = 1:650
    MeanFd(i) = norm(Fd(:,i))
end
mean1=mean(MeanFd);

MeanFdrag = zeros(650,10);
for i = 1:10
    for j = 1:650
        MeanFdrag(j,i) = norm(Fdrag(:,j,i))
    end
end
mean2 = mean(MeanFdrag, 1)

Fdragpick = Fdrag(:,:,40);
for i = 1:650
    MeanFd1(i) = norm(Fdragpick(:,i))
end
mean3=mean(MeanFd1);















%%

random_numbers = rand(1, 100);

% Scale and shift to the range 0.3 to 0.7
scaled_random_numbers = 0.3 + (0.7 - 0.3) * random_numbers;
random_numbers1 = rand(1, 100);

% Scale and shift to the range 0.3 to 0.7
scaled_random_numbers1 = 0.3 + (0.7 - 0.3) * random_numbers1;

StartLoc = [scaled_random_numbers; scaled_random_numbers1];
for k = 1:100
    % Extract the k-th set of coefficients

    % Initialize temporary storage for x(t) for the current set of coefficients
    x_temp = zeros(3, length(t));
    
    % Compute x(t) for each time point
    for i = 1:length(t)
        x_temp(:, i) = calculateXold(t(i), 0, 0, 0, t(end), W, StartLoc(1,k), StartLoc(2,k), 0);
    end

    % Store the result in the 4D array
    X1(:, :, k) = x_temp;
end

Particles_locations = zeros(3, 651, 100);
Velocity_fields = Particles_locations;
Fdrag = zeros(3,650,100);
Totalenergy1 = zeros(650,1000);
Totalenergy = zeros(100,1);
Realdrags = Totalenergy;
Gap = zeros(3,650,100);
Vs = Gap;
for i = 1:100
Particles_locations(:,:,i) = get_particle_cell_Location(X1(:,:,i),64);
Velocity_fields(:,:,i) = get_vectors_at_locations(uField, vField, wField,Particles_locations(:,:,i));
Gap(:,:,i) = diff(X1(:,:,i), 1, 2);
Vs(:,:,i) = Gap(:,:,i)/dt;
    for j = 1:650
        Fdrag(:,j,i)=(Velocity_fields(:,j,i)-Vs(:,j,i))*norm(Velocity_fields(:,j,i)-Vs(:,j,i));
        Totalenergy1(j,i) = dot(Fdrag(:,j,i), Gap(:,j,i));
    end

    Totalenergy(i) = sum(Totalenergy1(:,i));
    Realdrags(i)=Totalenergy(i)/(0.650);


end
bbb = Totalenergy/(0.650);
histogram(bbb/Realdrag)
title('Histogram of 100 strating point');
xlabel('Normalized COT');
ylabel('Frequency');
grid on;  % Add a grid for better readability





%%

% Define the time array and other parameters
dt = 5e-4;
t = 0:dt:650*dt;
W = 2;
StartLoc = [0.5, 0.5];

% Initial guess for the Fourier coefficients
initial_coeffs = zeros(12, 1);  % Example initial coefficients (4 coefficients for each of a, b, c)


% Call objective function directly to check for errors
totalEnergy = objective_function(optimized_coeffs, t, W, StartLoc, uField, vField, wField, dt);
disp('Initial total energy:');
disp(totalEnergy);

% Call objective function directly to check for errors
totalEnergy = objective_function(initial_coeffs, t, W, StartLoc, uField, vField, wField, dt);
disp('Initial total energy:');
disp(totalEnergy);




%%

% Define the time array and other parameters


% % Initial guess for the Fourier coefficients
% initial_coeffs = zeros(12, 1);  % Example initial coefficients (10 coefficients for each of a, b, c)
% 
% % Optimization options
% options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter');
% 
% % Perform optimization using fminunc
% result = fminunc(@(coeffs) objective_function(coeffs, t, W, StartLoc, uField, vField, wField, dt), initial_coeffs, options);
% 
% % Optimized coefficients
% optimized_coeffs = result;
% 
% disp('Optimized Fourier Coefficients:');
% disp(optimized_coeffs);

%%

% Define bounds for the coefficients
lb = -0.1 * ones(12, 1);  % Lower bounds
ub = 0.1*ones(12, 1);   % Upper bounds

% Optimization options
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter');

% Perform optimization using fmincon
result = fmincon(@(coeffs) objective_function(coeffs, t, W, StartLoc, uField, vField, wField, dt), ...
                 initial_coeffs, [], [], [], [], lb, ub, [], options);

% Optimized coefficients
optimized_coeffs = result;

disp('Optimized Fourier Coefficients:');
disp(optimized_coeffs);


%%

TKE = calculateTKE(uField, vField, wField);
disp(TKE);

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
    

% Define bounds for the coefficients
lb = -0.1 * ones(12, 1);  % Lower bounds
ub = 0.1*ones(12, 1);   % Upper bounds

% Optimization options
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter');

result = fmincon(@(coeffs) objective_function(coeffs, t, W, StartLoc, uFieldMultiplied, vFieldMultiplied, wFieldMultiplied, dt), ...
                 initial_coeffs, [], [], [], [], lb, ub, [], options);

% Optimized coefficients
optimized_coeffs = result;


TKE = calculateTKE(uFieldMultiplied, vFieldMultiplied, wFieldMultiplied);
disp(TKE);

%%

% Define bounds for the coefficients
lb = -0.1 * ones(12, 1);  % Lower bounds
ub = 0.1*ones(12, 1);   % Upper bounds

% Optimization options
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter');

% Perform optimization using fmincon
result = fmincon(@(coeffs) objective_function(coeffs, t, W, StartLoc, uFieldMultiplied, vFieldMultiplied, wFieldMultiplied, dt), ...
                 initial_coeffs, [], [], [], [], lb, ub, [], options);


% Optimized coefficients
optimized_coeffs = result;

disp('Optimized Fourier Coefficients:');
disp(optimized_coeffs);
% Call objective function directly to check for errors
totalEnergy = objective_function(optimized_coeffs, t, W, StartLoc, uFieldMultiplied, vFieldMultiplied, wFieldMultiplied, dt);
disp('Initial total energy:');
disp(totalEnergy);


%%

Tra = calculateX(t, optimized_coeffs,W, StartLoc);


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
plot3(Tra(1,:), Tra(2,:),Tra(3,:))



% Set specific axes limits
xlim([0 1]);  % X-axis limits
ylim([0 1]);  % Y-axis limits
zlim([0 1]);  % Z-axis limits

% Adjust the view angle for better 3D perception
view(3); % Default 3D view

hold off; % Release the hold on the current figure

Ploc = get_particle_cell_Location(Tra, 64);
Floc = get_vectors_at_locations(uField, vField, wField, Ploc);


% Extracting the components
X = Ploc(1, :);
Y = Ploc(2, :);
Z = Ploc(3, :);

U = Floc(1, :);
V = Floc(2, :);
W = Floc(3, :);

% Scale the vectors to make them longer
scaleFactor = 2; % Adjust this factor to make vectors longer or shorter
U = U * scaleFactor;
V = V * scaleFactor;
W = W * scaleFactor;

% Plotting the vectors in 3D
figure;
quiver3(X, Y, Z, U, V, W, 'LineWidth', 1.5);
hold on;
plot3(X, Y, Z, 'k-', 'LineWidth', 1); % Plot lines connecting the locations

xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Vector Plot');
grid on;
axis equal;


hold off;





