numFields = 5;     % You have fields1 to fields5
start_num = 350;
end_num = 630;
numb = end_num-start_num;
dt = 1e-3;
t = 0:dt:numb*dt;
W = 8;
p = 1;
tau_p = 1e-2;
Forcing = 1;
N = 13; % Number of parallel searches
% resultsU8ul5 = cell(1, N);
% fvalsU8ul5 = cell(1, N);
% 
% resultsTest1 = cell(1, N);
% fvalsTest1 = cell(1, N);
% 
% 
% resultsU8ul1high = cell(1, N);
% fvalsU8ul1high = cell(1, N);

%U8_3 t=150
%%
tic

G = [0.001,0.01,0.1,0.2,0.4,0.5,0.6,0.7,0.9,1.0,1.2,1.7,2];

Ghigh = [2.5,3,4,5,6,7,8,9,10,11,15,20,25];

% Initialize result storage
resultsU8_5 = cell(numFields, 1);
fvalsU8_5   = cell(numFields, 1);
resultsU8_10 = cell(numFields, 1);
fvalsU8_10   = cell(numFields, 1);
resultsU8_15 = cell(numFields, 1);
fvalsU8_15   = cell(numFields, 1);
resultsU8_20 = cell(numFields, 1);
fvalsU8_20   = cell(numFields, 1);



% Initialize result storage
resultsU8_5high = cell(numFields, 1);
fvalsU8_5high   = cell(numFields, 1);
resultsU8_10high = cell(numFields, 1);
fvalsU8_10high   = cell(numFields, 1);
resultsU8_15high = cell(numFields, 1);
fvalsU8_15high   = cell(numFields, 1);
resultsU8_20high = cell(numFields, 1);
fvalsU8_20high   = cell(numFields, 1);


for fieldIdx = 1:numFields
    % Define the file path dynamically
    file_path = sprintf('/Users/boyi/Simulations/hit/output_fields%dU=10', fieldIdx);

    % Read all field files
    numFiles = end_num - start_num + 1;
    variables_list = cell(numFiles, 1);
    
    for i = start_num:end_num
        filename = fullfile(file_path, sprintf('field_%05d', i));
        [~, ~, ~, ~, variables] = read_field_file(filename);
        variables_list{i - start_num + 1} = variables;
    end
    
    % Extract u, v, w fields
    uField = cell(numFiles, 1);
    vField = cell(numFiles, 1);
    wField = cell(numFiles, 1);
    
    for i = 1:numFiles
        data = variables_list{i};
        uField{i} = data(1);
        vField{i} = data(2);
        wField{i} = data(3);
    end
    
    % Calculate RMS velocity U
    U = calculateRMS(uField(1:100), vField(1:100), wField(1:100));
    
    % First optimization
    tempResults5 = cell(1, N);
    tempFval5 = cell(1, N);
    parfor i = 1:N
        [tempResults5{i}.coeffs, tempResults5{i}.W, tempResults5{i}.energy, tempFval5{i}] = ...
            optimization27(t, U, uField, vField, wField, dt, p, U, Forcing, 1000, 5, G(i), 2, 30, 1, 1, []);
    end
    resultsU8_5{fieldIdx} = tempResults5;
    fvalsU8_5{fieldIdx} = tempFval5;
    
    % Second optimization
    tempResults10 = cell(1,N);
    tempFval10 = cell(1,N);
    parfor i = 1:N
        [tempResults10{i}.coeffs, tempResults10{i}.W, tempResults10{i}.energy, tempFval10{i}] = ...
            optimization27(t, tempResults5{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000, 10, G(i), 2, 30, 1, 1, tempResults5{i}.coeffs);
    end
    resultsU8_10{fieldIdx} = tempResults10;
    fvalsU8_10{fieldIdx} = tempFval10;

    tempResults15 = cell(1, N);
tempFval15 = cell(1, N);
parfor i = 1:N
    [tempResults15{i}.coeffs, tempResults15{i}.W, tempResults15{i}.energy, tempFval15{i}] = ...
        optimization27(t, tempResults10{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000, 15, G(i), 2, 30, 1, 1, tempResults10{i}.coeffs);
end
resultsU8_15{fieldIdx} = tempResults15;
fvalsU8_15{fieldIdx} = tempFval15;

% Fourth optimization (step = 20)
tempResults20 = cell(1, N);
tempFval20 = cell(1, N);
parfor i = 1:N
    [tempResults20{i}.coeffs, tempResults20{i}.W, tempResults20{i}.energy, tempFval20{i}] = ...
        optimization27(t, tempResults15{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000, 20, G(i), 2, 30, 1, 1, tempResults15{i}.coeffs);
end
resultsU8_20{fieldIdx} = tempResults20;
fvalsU8_20{fieldIdx} = tempFval20;








% First optimization
tempResults5 = cell(1, N);
tempFval5 = cell(1, N);
parfor i = 1:N
    [tempResults5{i}.coeffs, tempResults5{i}.W, tempResults5{i}.energy, tempFval5{i}] = ...
        optimization27(t, U, uField, vField, wField, dt, p, U, Forcing, 1000, 5, Ghigh(i), 2, 300, 1, 1, []);
end
resultsU8_5high{fieldIdx} = tempResults5;
fvalsU8_5high{fieldIdx} = tempFval5;

% Second optimization
tempResults10 = cell(1,N);
tempFval10 = cell(1,N);
parfor i = 1:N
    [tempResults10{i}.coeffs, tempResults10{i}.W, tempResults10{i}.energy, tempFval10{i}] = ...
        optimization27(t, tempResults5{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000, 10, Ghigh(i), 2, 300, 1, 1, tempResults5{i}.coeffs);
end
resultsU8_10high{fieldIdx} = tempResults10;
fvalsU8_10high{fieldIdx} = tempFval10;

% Third optimization (step = 15)
tempResults15 = cell(1, N);
tempFval15 = cell(1, N);
parfor i = 1:N
    [tempResults15{i}.coeffs, tempResults15{i}.W, tempResults15{i}.energy, tempFval15{i}] = ...
        optimization27(t, tempResults10{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000, 15, Ghigh(i), 2, 300, 1, 1, tempResults10{i}.coeffs);
end
resultsU8_15high{fieldIdx} = tempResults15;
fvalsU8_15high{fieldIdx} = tempFval15;

% Fourth optimization (step = 20)
tempResults20 = cell(1, N);
tempFval20 = cell(1, N);
parfor i = 1:N
    [tempResults20{i}.coeffs, tempResults20{i}.W, tempResults20{i}.energy, tempFval20{i}] = ...
        optimization27(t, tempResults15{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000, 20, Ghigh(i), 2, 300, 1, 1, tempResults15{i}.coeffs);
end
resultsU8_20high{fieldIdx} = tempResults20;
fvalsU8_20high{fieldIdx} = tempFval20;

end


resultsall = cell(numFields,1);
fvalsall = resultsall;
Gall = [G, Ghigh];
for i = 1:numFields
    resultsall{i} = [resultsU8_20{i}, resultsU8_20high{i}];
    fvalsall{i} = [fvalsU8_20{i}, fvalsU8_20high{i}];% now each is 1x26
end




toc


%%



fvalue = [fvalsallU4_1; fvalsallU4_2;fvalsallU4_3; fvalsallU4_4];
%fvalue = [fvalsallU8_1; fvalsallU8_2; fvalsallU8_4];
N = numel(fvalue{1});  % Assuming each is a cell array of size N×1
M = numel(fvalue);
vals = zeros(M, N);  
for i = 1:N
   
    for fieldIdx = 1:M
        vals(fieldIdx, i) = fvalue{fieldIdx}{i};
    end
    
end
meanFvalsall = mean(vals, 1);
stdFvalsall = std(vals);

results = [resultsallU4_1; resultsallU4_2;resultsallU4_3; resultsallU4_4];
%results = [resultsallU8_1; resultsallU8_2; resultsallU8_4];

Wvals = zeros(M, N);  % 5 fields × 13 values
for i = 1:N
   
    for fieldIdx = 1:M
        Wvals(fieldIdx, i) = results{fieldIdx}{i}.W;  % extract 1×13 vector
    end
    
end
meanWall = mean(Wvals, 1);  % average across 5 rows → 1×13 vector

%%

% Preallocate E_min
E_min = zeros(size(Gall));
E_min_2 = zeros(size(Gall));

for i = 1:length(Gall)
        E_min(i) = EF_min_general(1, 0, 0,Gall(i));
        E_min_2(i) = EF_min_general(1.6, 0.6, 0.45,Gall(i));
end

% Plot
figure;
hold on;

defaultColor = [0 0 0]; % MATLAB default blue
plot(Gall, E_min / 1.611, '-', 'LineWidth', 2, 'Color', defaultColor);
plot(Gall, E_min_2 / 1.611, '-', 'LineWidth', 2, 'Color', [1 0 0]);


% Preallocate arrays
fvalsp = zeros(size(Gall));



% Convert cell to array
for i = 1:length(Gall)
    fvalsp(i) = meanFvalsall(i);
end

% Define consistent colors (can adjust as you like)
color1 = [0 0.4470 0.7410]; % MATLAB blue
color2 = [0.8500 0.3250 0.0980]; % MATLAB reddish

% Plot fvals1 on G and Ghigh, same marker & color
errorbar(Gall, fvalsp,stdFvalsall, '--s', 'LineWidth', 2, 'Color', color1, 'MarkerSize', 8);


% Axes and formatting
set(gca, 'XScale', 'log');
xlabel('G');
ylabel('E_{tailwind}');
title('Energy vs G (B = 1)');
grid on;

% Legend with larger font
legend( 'E_{TW} tailwind=1 updraft=0, crosswind = 0', 'E_{TW} tailwind=1.6 updraft=0.6, crosswind = 0.5', 'OptE trial 1', ...
       'Location', 'best', 'FontSize', 14);