x = [10 5 3 2];
x = [3 5 10];
col = 4;

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

figure
hold on

colors = lines(4);

legendNames = { ...
    'G = 0.2', ...
    'G = 1', ...
    'G = 25', ...
};

tn = [3,9,25];

for i = 1:3
    scatter(updraft_avgstotal(:,tn(i)), ...
            ywind_avgstotal(:,tn(i)), ...
            120, ...            % bigger marker size
            colors(i,:), ...
            'filled', ...
            'DisplayName', legendNames{i});
end

theta = linspace(0,2*pi,300);
plot(1.2*cos(theta), 1.2*sin(theta), 'k--', 'LineWidth',1.5, ...
     'DisplayName','r = 1.2')

xlabel('Updraft')
ylabel('Tail')
axis equal
grid on
xlim([0 2])
ylim([0 2])

legend('Location','best')
hold off




