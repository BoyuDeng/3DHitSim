parfor i = 1:N
    [results2Z4_5{i}.coeffs, results2Z4_5{i}.W, results2Z4_5{i}.energy, fvals2Z4_5{i}] = ...
        optimization27(t, U, uField, vField, wField, dt, p, U, Forcing, 1000,5, G(i),0, 20,1,1,0);
end 

parfor i = 1:N
    [results2Z4_10{i}.coeffs, results2Z4_10{i}.W, results2Z4_10{i}.energy, fvals2Z4_10{i}] = ...
        optimization27(t, results2Z4_5{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000,10, G(i),0, 20,1,1,results2Z4_5{i}.coeffs);
end

parfor i = 1:N
    [results2Z4_15{i}.coeffs, results2Z4_15{i}.W, results2Z4_15{i}.energy, fvals2Z4_15{i}] = ...
        optimization27(t, results2Z4_10{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000,15, G(i),0, 20,1,1,results2Z4_10{i}.coeffs);
end

parfor i = 1:N
    [results2Z4_20{i}.coeffs, results2Z4_20{i}.W, results2Z4_20{i}.energy, fvals2Z4_20{i}] = ...
        optimization27(t, results2Z4_15{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000,20, G(i),0, 20,1,1,results2Z4_15{i}.coeffs);
end

parfor i = 1:N
    [results2Z4_25{i}.coeffs, results2Z4_25{i}.W, results2Z4_25{i}.energy, fvals2Z4_25{i}] = ...
        optimization27(t, results2Z4_20{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000,25, G(i),0, 20,1,1,results2Z4_20{i}.coeffs);
end
%%

%%
% Given values
%G = [ 0.1; 0.2; 0.4; 0.5; 0.6; 0.7; 0.9; 1.0; 1.2; 1.7; 2; 2.5];
G = [0.001;0.01;0.1;0.2;0.4;0.5;0.6;0.7;0.9;1.0;1.2;1.7;2;2.5];
B = 1;
c = 0.5;

% Preallocate E_min
E_min_06 = zeros(size(G));
E_min_09 = zeros(size(G));
u_crit = zeros(size(G));

% Define the expression for the critical point (u)
for i = 1:length(G)
% Define the function f(u)
u_crit(i) = (-1 + sqrt(1 - 4*(4*c*G(i) - 2*c^2 - 2*G(i)^2 - 2)))/2;
E_min_06(i) = (G(i)/u_crit(i)) * (1 + (1 + u_crit(i)^2 - 2*u_crit(i) + c^2)/G(i)^2 - (2*c/G(i)))^(3/4);

u_crit(i) = (-1 + sqrt(1 - 4*(4*c*Ghigh(i) - 2*c^2 - 2*Ghigh(i)^2 - 2)))/2;
E_min_high(i) = (Ghigh(i)/u_crit(i)) * (1 + (1 + u_crit(i)^2 - 2*u_crit(i) + c^2)/Ghigh(i)^2 - (2*c/Ghigh(i)))^(3/4);

% 
% E_min_09(i) = EF_min_07(G(i));
% E_min_high(i) = EF_min_07(Ghigh(i));
end

% Plot
figure;
hold on;

% Plot normalized E_min for c = 0.6 (no markers)
%plot(G, E_min_06 / 1.611, '-', 'LineWidth', 2);

% Plot normalized E_min for c = 0.3 (no markers)
plot(G, E_min_06 / 1.611, '-', 'LineWidth', 2);
plot(Ghigh, E_min_high / 1.611, '-', 'LineWidth', 2);

% Preallocate arrays
fvals2Z4 = zeros(size(G));
fvals2U4 = zeros(size(G));
fvalshZ4 = zeros(size(G));

% Convert cell to array
for i = 1:length(G)
    fvals2Z4(i) = fvals2Z4_10{i};
    %fvals2U4(i) = fvals2U4_25{i};
    fvalshZ4(i) = fvalshighZ4_10{i};
end

% Plot all three OptE datasets with markers
plot(G, fvals2Z4, '--s', 'LineWidth', 2);
%plot(G, fvals2U4, '--d', 'LineWidth', 2);
plot(Ghigh, fvalshZ4, '--^', 'LineWidth', 2);

% Axes and formatting
set(gca, 'XScale', 'log');
xlabel('G');
ylabel('E_{tailwind}');
title('Energy vs G (B = 1)');
grid on;

% Legend with larger font
legend( 'E_{TW} wy^2=1 wz^2 = 1', 'E_{TW} wy^2=1 wz^2 = 1', 'OptE U=4', 'OptE U=4', ...
       'Location', 'best', 'FontSize', 14);

%%

% G values
G = [0.001, 0.01, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.9, 1.0, 1.2, 1.7, 2, 2.5];

% Compute minimum E using the three functions
E07 = EF_min_07(G);
E08 = EF_min_08(G);
E09 = EF_min_09(G);

% Plot the results
figure;
plot(G, E09/1.611, '-o', 'LineWidth', 1.5); hold on;
plot(G, E08/1.611, '-s', 'LineWidth', 1.5);
plot(G, E07/1.611, '-^', 'LineWidth', 1.5);
plot(G, E_min_06/1.611, '-r', 'LineWidth', 1.5);
grid on;

% Labeling
set(gca, 'XScale', 'log');
xlabel('G');
ylabel('Minimum E');
title('Comparison of E\_min for Different y^2 and c^2');
legend('EF\_min\_09 (y^2=0.9, c^2=0.1)', ...
       'EF\_min\_08 (y^2=0.8, c^2=0.2)', ...
       'EF\_min\_07 (y^2=0.7, c^2=0.3)', ...
       'EF\_min\_11 (y^2=1, c^2=1)', ...
       'Location', 'northwest', ...
       'FontSize', 18);  % Adjust font size as needed

