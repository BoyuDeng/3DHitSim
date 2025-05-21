colors = {'b', 'r', 'g', 'm', 'c', 'k'}; % Extend this if needed
markers = {'o', 's', 'd', '^', 'v', 'p'}; % Extend this if needed

hold on; % Hold the plot to overlay multiple points

for i = 1:length(G)
    % Plot each fvals dataset with a different color and marker
    plot(G(i), fvals1U8_5{i}, markers{1}, 'MarkerSize', 6, 'MarkerFaceColor', colors{1});
    plot(G(i), fvals1U8_10{i}, markers{2}, 'MarkerSize', 6, 'MarkerFaceColor', colors{2});
    plot(G(i), fvals1U8_15{i}, markers{3}, 'MarkerSize', 6, 'MarkerFaceColor', colors{3});
    plot(G(i), fvals1U8_20{i}, markers{4}, 'MarkerSize', 6, 'MarkerFaceColor', colors{4});
    plot(G(i), fvals1U8_25{i}, markers{5}, 'MarkerSize', 6, 'MarkerFaceColor', colors{5});
end

set(gca, 'XScale', 'log'); % Set X-axis to log scale
set(gca, 'YScale', 'log'); % Set X-axis to log scale
% Adding labels and title
xlabel('G');
ylabel('COT');
title('Plot of COT vs G t = 16.5 (t_d U/L)');

grid on;

% Adding legend
legend({'mode5', 'mode10', 'mode15','mode20','mode25'}, 'Location', 'best');

hold off; % Release the plot hold

%%
colors = {'b', 'r', 'g', 'm', 'c', 'k'}; % Extend this if needed
markers = {'o', 's', 'd', '^', 'v', 'p'}; % Extend this if needed

hold on; % Hold the plot to overlay multiple points

for i = 1:length(G)
    % Plot each fvals dataset with a different color and marker
    plot(G(i), fvals1U4_5{i}, markers{1}, 'MarkerSize', 6, 'MarkerFaceColor', colors{1});
    plot(G(i), fvals1U4_10{i}, markers{2}, 'MarkerSize', 6, 'MarkerFaceColor', colors{2});
    plot(G(i), fvals1U4_15{i}, markers{3}, 'MarkerSize', 6, 'MarkerFaceColor', colors{3});
    plot(G(i), fvals1U4_20{i}, markers{4}, 'MarkerSize', 6, 'MarkerFaceColor', colors{4});
    plot(G(i), fvals1U4_25{i}, markers{5}, 'MarkerSize', 6, 'MarkerFaceColor', colors{5});
end

set(gca, 'XScale', 'log'); % Set X-axis to log scale
set(gca, 'YScale', 'log'); % Set Y-axis to log scale

% Adding labels and title
xlabel('G');
ylabel('COT');
title('Plot of COT vs G U_ud=4, t = 16.5');

grid on;

% Adding legend
legend({'mode5', 'mode10', 'mode15','mode20','mode25'}, 'Location', 'best');

hold off; % Release the plot hold

%%

% Calculate improvements (if not done yet)
improvement_5to10 = zeros(1, length(G));
improvement_10to15 = zeros(1, length(G));
improvement_15to20 = zeros(1, length(G));
improvement_20to25 = zeros(1, length(G));

for i = 1:length(G)
    improvement_5to10(i)  = (fvals1U4_5{i}  - fvals1U4_10{i})  / fvals1U4_5{i};
    improvement_10to15(i) = (fvals1U4_10{i} - fvals1U4_15{i}) / fvals1U4_10{i};
    improvement_15to20(i) = (fvals1U4_15{i} - fvals1U4_20{i}) / fvals1U4_15{i};
    improvement_20to25(i) = (fvals1U4_20{i} - fvals1U4_25{i}) / fvals1U4_20{i};
end

% Plotting
figure;
hold on;

plot(G, improvement_5to10, '-o', 'LineWidth', 1.5);
plot(G, improvement_10to15, '-s', 'LineWidth', 1.5);
plot(G, improvement_15to20, '-^', 'LineWidth', 1.5);
plot(G, improvement_20to25, '-d', 'LineWidth', 1.5);

set(gca, 'XScale', 'log'); % Optional: log scale for G
xlabel('G');
ylabel('Relative Improvement');
title('Relative Improvement of COT between Mode Increments');
legend({'5→10', '10→15', '15→20', '20→25'}, 'Location', 'best');
grid on;
hold off;


%%
% Given values
G = [0.001; 0.01; 0.1; 0.2; 0.4; 0.5; 0.6; 0.7; 0.9; 1.0; 1.2; 1.7; 2; 2.5;10;30;200];
B = 1;

% Preallocate E_min
E_min = zeros(size(G));

for i = 1:length(G)
    Gi = G(i);
    
    % Compute U*
    U_star = 0.5 * (sqrt(8 * Gi^2 * B^2 + 9) - 1);
    
    % Compute the energy at U*
    numerator = ((1 - U_star)^2 + Gi^2 * B^2)^(3/4);
    denominator = U_star * sqrt(Gi);
    E_min(i) = numerator / denominator;
end

% Plot
figure;
plot(G, E_min/1.611, '-o', 'LineWidth', 2);
set(gca, 'XScale', 'log'); % Optional: log scale for G
xlabel('G');
ylabel('Minimum E_{TW}');
title('Minimum Energy vs G (B = 1)');
grid on;


%%
% Given values
%G = [ 0.1; 0.2; 0.4; 0.5; 0.6; 0.7; 0.9; 1.0; 1.2; 1.7; 2; 2.5];
G = [0.001;0.01;0.1;0.2;0.4;0.5;0.6;0.7;0.9;1.0;1.2;1.7;2;2.5];
B = 1;
c_1 = 0.6;
c_2 = 0;

% Preallocate E_min
E_min_06 = zeros(size(G));

for i = 1:length(G)
    Gi = G(i);
    
    % Compute U* based on modified formula with c^2 = 1
    U_star = 0.5 * (sqrt(8 * Gi^2 * B^2 + 8 * c_1^2 + 9) - 1);
    
    % Compute energy at U*
    A = (1 - U_star)^2 + Gi^2 * B^2 + c_1^2;
    numerator = A^(3/4);
    denominator = U_star * sqrt(Gi);
    E_min_06(i) = numerator / denominator;
end

E_min_03 = zeros(size(G));

for i = 1:length(G)
    Gi = G(i);
    
    % Compute U* based on modified formula with c^2 = 1
    U_star = 0.5 * (sqrt(8 * Gi^2 * B^2 + 8 * c_2^2 + 9) - 1);
    
    % Compute energy at U*
    A = (1 - U_star)^2 + Gi^2 * B^2 + c_2^2;
    numerator = A^(3/4);
    denominator = U_star * sqrt(Gi);
    E_min_03(i) = numerator / denominator;
end

% Plot
figure;
hold on;

% Plot normalized E_min for c = 0.6 (no markers)
%plot(G, E_min_06 / 1.611, '-', 'LineWidth', 2);

% Plot normalized E_min for c = 0.3 (no markers)
plot(G, E_min_03 / 1.611, '-', 'LineWidth', 2);

% Preallocate arrays
fvals2Z4 = zeros(size(G));
fvals2U4 = zeros(size(G));
fvals1U8 = zeros(size(G));

% Convert cell to array
for i = 1:length(G)
    fvals2Z4(i) = fvals2Z4_10{i};
    fvals2U4(i) = fvals2U4_25{i};
    fvals1U8(i) = fvals1U8_25{i};
end

% Plot all three OptE datasets with markers
plot(G, fvals2Z4, '--s', 'LineWidth', 2);
plot(G, fvals2U4, '--d', 'LineWidth', 2);
plot(G, fvals1U8, '--^', 'LineWidth', 2);

% Axes and formatting
set(gca, 'XScale', 'log');
xlabel('G');
ylabel('E_{tailwind}');
title('Energy vs G (B = 1)');
grid on;

% Legend with larger font
legend( 'E_{TW} 2D flow', 'OptE U=4 no Potential Change', 'OptE U=4', 'OptE U=8', ...
       'Location', 'best', 'FontSize', 14);




%%

colors = {'b', 'r', 'g', 'm', 'c', 'k'}; % Extend this if needed
markers = {'o', 's', 'd', '^', 'v', 'p'}; % Extend this if needed

hold on; % Hold the plot to overlay multiple points

for i = 1:length(G)
    % Plot each fvals dataset with a different color and marker
    plot(G(i), fvals2U4_5{i}, markers{1}, 'MarkerSize', 6, 'MarkerFaceColor', colors{1});
    plot(G(i), fvals2U4_10{i}, markers{2}, 'MarkerSize', 6, 'MarkerFaceColor', colors{2});
    plot(G(i), fvals2U4_15{i}, markers{3}, 'MarkerSize', 6, 'MarkerFaceColor', colors{3});
    plot(G(i), fvals2U4_20{i}, markers{4}, 'MarkerSize', 6, 'MarkerFaceColor', colors{4});
    plot(G(i), fvals2U4_25{i}, markers{5}, 'MarkerSize', 6, 'MarkerFaceColor', colors{5});
end

set(gca, 'XScale', 'log'); % Set X-axis to log scale
set(gca, 'YScale', 'log'); % Set Y-axis to log scale

% Adding labels and title
xlabel('G');
ylabel('COT');
title('Plot of COT vs G U_ud=4, t = 16.5');

grid on;

% Adding legend
legend({'mode5', 'mode10', 'mode15','mode20','mode25'}, 'Location', 'best');

hold off; % Release the plot hold

