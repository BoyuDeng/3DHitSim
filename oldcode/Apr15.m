parfor i = 1:N
    [results2U4_5{i}.coeffs, results2U4_5{i}.W, results2U4_5{i}.energy, fvals2U4_5{i}] = ...
        optimization27(t, U, uField, vField, wField, dt, p, U, Forcing, 1000,5, G(i),0, 20,1,1,0);
end 

parfor i = 1:N
    [results2U4_10{i}.coeffs, results2U4_10{i}.W, results2U4_10{i}.energy, fvals2U4_10{i}] = ...
        optimization27(t, results2U4_5{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000,10, G(i),0, 20,1,1,results2U4_5{i}.coeffs);
end

parfor i = 1:N
    [results2U4_15{i}.coeffs, results2U4_15{i}.W, results2U4_15{i}.energy, fvals2U4_15{i}] = ...
        optimization27(t, results2U4_10{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000,15, G(i),0, 20,1,1,results2U4_10{i}.coeffs);
end

parfor i = 1:N
    [results2U4_20{i}.coeffs, results2U4_20{i}.W, results2U4_20{i}.energy, fvals2U4_20{i}] = ...
        optimization27(t, results2U4_15{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000,20, G(i),0, 20,1,1,results2U4_15{i}.coeffs);
end

parfor i = 1:N
    [results2U4_25{i}.coeffs, results2U4_25{i}.W, results2U4_25{i}.energy, fvals2U4_25{i}] = ...
        optimization27(t, results2U4_20{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000,25, G(i),0, 20,1,1,results2U4_20{i}.coeffs);
end



%%
COT = zeros(1, length(G));

ali = zeros(1, length(G));
for i = 1:14
[COT(i),~, ~, ali(i)] = COT14(results1U4_25{i}.coeffs, t, results1U4_25{i}.W, uField, vField, wField, dt, p,U, 1, G(i));
%ali(i) = sum(vecnorm(Fdrag(:,:,i)))/700;
end

%%
% Given G values
%G = [0.1; 0.2; 0.4; 0.5; 0.6; 0.7; 0.9; 1.0; 1.2; 1.7; 2; 2.5];
B = 1;
c = 0;

% Preallocate energy array
E_min = zeros(size(G));

for i = 1:length(G)
    Gi = G(i);
    
    % Compute U* for fixed c
    U_star = 0.5 * (sqrt(8 * Gi^2 * B^2 + 8 * c^2 + 9) - 1);
    
    % Compute energy
    A = (1 - U_star)^2 + Gi^2 * B^2 + c^2;
    E_min(i) = A^(3/4) / (U_star * sqrt(Gi));
end

% Plot
figure;
semilogx(G, E_min/1.611, '-o', 'LineWidth', 2);
xlabel('G');
ylabel('Minimum E_{TW,c}');
title('Minimum Energy vs G (c = 1, B = 1)');
grid on;
