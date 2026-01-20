x = [10 5 3 2];
col = 2;

figure
plot(x, updraft_avgstotal(:,col), '.', 'MarkerSize', 20)
hold on
plot(x, crosswind_avgstotal(:,col), '.', 'MarkerSize', 20)
plot(x, ywind_avgstotal(:,col), '.', 'MarkerSize', 20)
hold off

xlabel('X values')
ylabel('Data')
legend('Updraft', 'Crosswind', 'Y-wind')
grid on

figure
scatter(updraft_avgstotal(:), ywind_avgstotal(:), 'filled')
hold on

theta = linspace(0, 2*pi, 300);
r = 1.2;
plot(r*cos(theta), r*sin(theta), 'k--', 'LineWidth', 1.5)

xlabel('Updraft')
ylabel('Tail')
grid on
axis equal
xlim([0 2])
ylim([0 2])

hold off



