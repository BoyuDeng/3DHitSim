%%
colors = {'b', 'r', 'g', 'm', 'c', 'k'}; % Extend this if needed
markers = {'o', 's', 'd', '^', 'v', 'p'}; % Extend this if needed

hold on; % Hold the plot to overlay multiple points

for i = 1:length(G)
    % Plot each fvals dataset with a different color and marker
    plot(G(i), fvalsU8_35_m5{i}, markers{1}, 'MarkerSize', 6, 'MarkerFaceColor', colors{1});
    plot(G(i), fvalsU8_35_m10{i}, markers{2}, 'MarkerSize', 6, 'MarkerFaceColor', colors{2});
    plot(G(i), fvalsU8_35_m15{i}, markers{3}, 'MarkerSize', 6, 'MarkerFaceColor', colors{3});
    plot(G(i), fvalsU8_35_m20{i}, markers{4}, 'MarkerSize', 6, 'MarkerFaceColor', colors{4});
    plot(G(i), fvalsU8_35_m25{i}, markers{5}, 'MarkerSize', 6, 'MarkerFaceColor', colors{5});
end

set(gca, 'XScale', 'log'); % Set X-axis to log scale
set(gca, 'YScale', 'log'); % Set X-axis to log scale
% Adding labels and title
xlabel('G');
ylabel('COT');
title('Plot of COT vs G t = 18.6 (t_d U/L)');

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
    plot(G(i), fvalsU8_m30{i}, markers{1}, 'MarkerSize', 6, 'MarkerFaceColor', colors{1});
    plot(G(i), fvalsU8_35_m25{i}, markers{2}, 'MarkerSize', 6, 'MarkerFaceColor', colors{2});
    plot(G(i), fvalsU4_m25{i}, markers{3}, 'MarkerSize', 6, 'MarkerFaceColor', colors{3});
    plot(Ghigh(i), fvalsU8_35_m15high{i}, markers{4}, 'MarkerSize', 6, 'MarkerFaceColor', colors{4});
    %plot(G(i), fvalsU8_35_m25{i}, markers{5}, 'MarkerSize', 6, 'MarkerFaceColor', colors{5});
end

set(gca, 'XScale', 'log'); % Set X-axis to log scale
set(gca, 'YScale', 'log'); % Set X-axis to log scale
% Adding labels and title
xlabel('G');
ylabel('COT');
title('Plot of COT vs G');

grid on;

% Adding legend
legend({'t=37.3 m30', 't=18.6 m25', 't=16.5 m25','t=18.6 m15'}, 'Location', 'best');

hold off; % Release the plot hold

%%
for i = 1:length(G)
    X = X14unlim(t, results1U8_20{i}.coeffs, results1U8_20{i}.W);
    Zdiffh(i) = X(3, end)-X(3, 1);
end


colors = {'b', 'r', 'g', 'm', 'c', 'k'}; % Extend this if needed
markers = {'o', 's', 'd', '^', 'v', 'p'}; % Extend this if needed

hold on; % Hold the plot to overlay multiple points

for i = 1:length(G)
    % Plot each fvals dataset with a different color and marker
    plot(G(i), Zdiffh1U8(i), markers{1}, 'MarkerSize', 10, 'MarkerFaceColor', colors{1});
    plot(G(i), Zdiffh1U4(i), markers{2}, 'MarkerSize', 10, 'MarkerFaceColor', colors{2});
    plot(G(i), Zdiffh2U4(i), markers{3}, 'MarkerSize', 10, 'MarkerFaceColor', colors{3});

end

set(gca, 'XScale', 'log'); % Set X-axis to log scale
%set(gca, 'YScale', 'log'); % Set X-axis to log scale
% Adding labels and title
xlabel('G');
ylabel('Change in Z');
title('Plot of diff Z vs G');

grid on;

% Adding legend
legend({'flow U=8', 'flow1 U=4', 'flow2 U=4'}, 'Location', 'best');

hold off; % Release the plot hold