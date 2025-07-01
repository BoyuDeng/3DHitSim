parfor i = 1:N
    [results2shZ8_5{i}.coeffs, results2shZ8_5{i}.W, results2shZ8_5{i}.energy, fvalsh2sZ8_5{i}] = ...
        optimization27(t, U, uField, vField, wField, dt, p, U, Forcing, 1000,5, Ghigh(i),2,300,1,1,0);
end 

parfor i = 1:N
    [results2shZ8_10{i}.coeffs, results2shZ8_10{i}.W, results2shZ8_10{i}.energy, fvals2shZ8_10{i}] = ...
        optimization27(t, resultsh2sZ8_5{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000,10, Ghigh(i),2,300,1,1,resultsh2sZ8_5{i}.coeffs);
end
%%
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
Wstar = zeros(3,13);
Ustar = zeros(3,13);
G = [0.001;0.01;0.1;0.2;0.4;0.5;0.6;0.7;0.9;1.0;1.2;1.7;2];
Ghigh = [2.5,3,4,5,6,7,8,9,10,11,15,20,25];
B = 1;
c = 0.5;

% Preallocate E_min
E_min_06 = zeros(size(G));
E_min_09 = zeros(size(G));
E_min_high= zeros(size(G));
u_crit = zeros(size(G));
for i = 1:13
coeffs = results2Z4_10{i}.coeffs;
W = results2Z4_10{i}.W;


        X = X14unlim(t, coeffs, W);

        Velocity_fields = get_vel(uField,vField,wField,64,X);
        
        Vs = V14unlim(t, coeffs, W);

        Utemp = abs(Vs);
        Wtemp = (Velocity_fields);
        Wstar(:,i) = mean(Wtemp,2)./U;
        Ustar(:,i) = mean(Utemp,2)./U;





        E_min_09(i) = EF_min_general(1, 1, 0.5,G(i));
        E_min_high(i) = EF_min_general(1, 1, 0.5,Ghigh(i));
        %[~,Uc(i)] = EF_min_general(Wstar(2,i), Wstar(3,i),G(i));
end

% Plot
figure;
hold on;

% Plot normalized E_min for c = 0.6 (no markers)
%plot(G, E_min_06 / 1.611, '-', 'LineWidth', 2);

defaultColor = [0 0 0]; % MATLAB default blue
plot(G, E_min_09 / 1.611, '-', 'LineWidth', 2, 'Color', defaultColor);
plot(Ghigh, E_min_high / 1.611, '-', 'LineWidth', 2, 'Color', defaultColor);


% Preallocate arrays
fvals1 = zeros(size(G));
fvals2 = zeros(size(G));
fvals1h = zeros(size(G));
fvals2h = zeros(size(G));

% Convert cell to array
for i = 1:length(G)
    fvals1(i) = fvals2sZ8_10{i};
    fvals2(i) = fvals4Z4_10{i};
    fvals1h(i) = fvalsh2sZ8_5{i};
    fvals2h(i) = fvals4hZ4_10{i};
end

% Define consistent colors (can adjust as you like)
color1 = [0 0.4470 0.7410]; % MATLAB blue
color2 = [0.8500 0.3250 0.0980]; % MATLAB reddish

% Plot fvals1 on G and Ghigh, same marker & color
plot(G, fvals1, '--s', 'LineWidth', 2, 'Color', color1, 'MarkerSize', 8);
plot(Ghigh, fvals1h, '--s', 'LineWidth', 2, 'Color', color1, 'MarkerSize', 8);

% Plot fvals2 on G and Ghigh, same marker & color
plot(G, fvals2, '--d', 'LineWidth', 2, 'Color', color2, 'MarkerSize', 8);
plot(Ghigh, fvals2h, '--d', 'LineWidth', 2, 'Color', color2, 'MarkerSize', 8);


% Axes and formatting
set(gca, 'XScale', 'log');
xlabel('G');
ylabel('E_{tailwind}');
title('Energy vs G (B = 1)');
grid on;

% Legend with larger font
legend( 'E_{TW} y^2=1 c^2=1, s^2 = 0.25', 'E_{TW} y^2=1, z^2=1, s^2 = 0.25', 'OptE trial 1', 'OptE trial 1', 'OptE trial 2','OptE trial 2',...
       'Location', 'best', 'FontSize', 14);
        
%%

% Preallocate storage
numCases = 14;
D_allhigh = [];  % Will hold all Dterm columns

for i = 1:numCases
    [~, ~, ~, ~, Dterm] = COT14(resultshighZ4_10{i}.coeffs, t, resultshighZ4_10{i}.W, ...
                                uField, vField, wField, ...
                                dt, p, U, Forcing, Ghigh(i));
    
    % Collect Dterm columns
    D_allhigh(:, i) = mean(Dterm,2);
end


