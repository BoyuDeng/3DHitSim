zdot_avgstotal=[];
zupdraft_avgstotal=[];
zcrosswind_avgstotal=[];
zywind_avgstotal = [];
cotstraightall = [];
%%
%fvalue = [fvalsallU4_1; fvalsallU4_2;fvalsallU4_3; fvalsallU4_4];
%fvalue = [fvalsallU8_1; fvalsallU8_2; fvalsallU8_4; fvalsallU8_5];
%fvalue = [fvalsallU10_1;fvalsallU10_2;fvalsallU10_3];
%fvalue = [fvalsallU4_1;fvalsallU4_2;fvalsallU4_3;fvalsallU4_4];
%results = [resultsallU4_1; resultsallU4_2;resultsallU4_3; resultsallU4_4];
%results = [resultsallU8_1; resultsallU8_2; resultsallU8_4; resultsallU8_5];
%results = [resultsallU10_1; resultsallU10_2; resultsallU10_3];
%results = [resultsallU4_1;resultsallU4_2;resultsallU4_3;resultsallU4_4];
numFields = length(results); 

numSamples = length(results{1});    % e.g., 13 or so

dot_avgs = zeros(numFields, numSamples);  % Store dot_avg for each sample in each field
updraft_avgs = zeros(numFields, numSamples);
Vrms = zeros(numFields, numSamples);
 crosswind_avgs = zeros(numFields, numSamples);
 ywind_avgs = zeros(numFields, numSamples);
 cotstraight = zeros(numFields, numSamples);



for n = 1:numFields

    % 
    % fn = mod(n, 5);
    % if fn == 0
    %    fn = 5;
    % end



    % Define the file path dynamically
    file_path = sprintf('/Users/boyi/Simulations/hit/output_fields%dU=8', 100);

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
    for i = 1:numSamples
        tempresults = results{n}{i};
        %coeffs = -1 + 2 * rand(1, 5);
        coeffs = [0 0 0 0 0];

        Wset = linspace(U, 2*U, 26);
        W= Wset(i);
        
        cotstraight(n,i) = COT14(coeffs, t, tempresults.W, uField,vField,wField, dt, 1, U, 1, Gall(i));
        

        X = X14unlim(t, coeffs, W);

        Velocity_fields = get_vel(uField, vField, wField, 64, X) ./ U;

        Vs = V14unlim(t, coeffs, W) ./ U;
        norm = sqrt(sum(Vs.^2,1));
        Vunit = Vs./ norm;


        dot_all = sum(Velocity_fields .* Vunit,1);  % 1×T
        dot_avg = mean(dot_all);
        
        Velsq = Velocity_fields.^2;
        Velmean = mean(Velsq, 2);
        Vrms(n,i) = sqrt(sum(Velmean)/3);
        
        
        updraft_all = Velocity_fields(3,:);
        ywind_all  = Velocity_fields(2,:);
        crosswind_all = abs(Velocity_fields(1,:));
        updraft_avg=mean(updraft_all);
         ywind_avg  = mean(ywind_all); 
        crosswind_avg = mean(crosswind_all);
        dot_avgs(n, i) = dot_avg;
        updraft_avgs(n,i) = updraft_avg;
        ywind_avgs(n,i)     = ywind_avg;
        crosswind_avgs(n,i) = crosswind_avg;
    end
end

zdot_avgstotal = [zdot_avgstotal; dot_avgs];
zupdraft_avgstotal = [zupdraft_avgstotal; updraft_avgs];
zcrosswind_avgstotal = [zcrosswind_avgstotal; crosswind_avgs];
zywind_avgstotal   = [zywind_avgstotal; ywind_avgs]; 
cotstraightall = [cotstraightall; cotstraight];


%%
% Means
mean_dot_per_field = mean(zdot_avgstotal, 1);
mean_updraft = mean(zupdraft_avgstotal, 1);
mean_crosswind = mean(zcrosswind_avgstotal, 1);
mean_ywind = mean(zywind_avgstotal,  1);

% Standard Errors
std_dot_per_field = std(zdot_avgstotal, 0, 1) / sqrt(size(zdot_avgstotal, 1));
std_updraft = std(zupdraft_avgstotal, 0, 1) / sqrt(size(zupdraft_avgstotal, 1));
std_crosswind = std(zcrosswind_avgstotal, 0, 1) / sqrt(size(zcrosswind_avgstotal, 1));
std_ywind = std(zywind_avgstotal, 0, 1) / sqrt(size(zywind_avgstotal,1));
Ux = linspace(1,2,26);
% Plot 1: <u · w> (Dot Product)
figure;
errorbar(Ux, mean_dot_per_field, std_dot_per_field, 'o-', 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
xlabel('U(Normilized Particle Velocity)');
ylabel('<w · u_{unit}>');
title('Dot Product of Flow and Trajectory');
grid on;
print(gcf, 'dot_product_plot.pdf', '-dpdf', '-bestfit');
% Plot 2: Updraft
figure;
errorbar(Ux, mean_updraft, std_updraft, 's-', 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
xlabel('U(Normilized Particle Velocity)');
ylabel('<W_{z}>');
title('Mean Updraft');
grid on;

% Plot 3: Crosswind
figure;
errorbar(Ux, mean_crosswind, std_crosswind, 'd-', 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
xlabel('U(Normilized Particle Velocity)');
ylabel('<|W_{x}|>');
title('Mean Crosswind');
grid on;

% Plot 4: ywind
figure;
errorbar(Ux, mean_ywind, std_ywind, '^-', 'LineWidth', 1.5);
set(gca,'XScale','log');
xlabel('U(Normilized Particle Velocity)');
ylabel('<W_{y}>');
title('Mean Y-wind');
grid on;

%%
% Suppose your cell array is named 'C'
total_sum = 0;
total_elements = 0;

for i = 1:350
    A = vField{i}; 
    A=A{1};% Extract the 63x63x63 array
    total_sum = total_sum + sum(A(:));  % Sum all elements
    total_elements = total_elements + numel(A);  % Count total elements
end

meancell = total_sum / total_elements;

%%

%% --- 1. Collect all values ------------------------------------------------
allVals = zeros(350 * 63^3, 1);
idx = 1;

for k = 1:350
    A  = vField{k}{1};             % 63×63×63 numeric array
    nA = numel(A);
    allVals(idx:idx+nA-1) = A(:) ./ U;  % normalize by U
    idx = idx + nA;
end
allVals = allVals(1:idx-1);

%% --- 2. Histogram as PDF --------------------------------------------------
[counts, edges] = histcounts(allVals, 'BinMethod', 'auto', ...
                             'Normalization', 'pdf');
centers = edges(1:end-1) + diff(edges)/2;

% Plot PDF
figure;
stairs(centers, counts, 'color', color1, 'LineWidth', 1.5);   % black outline
hold on;
grid on;

xlabel('w_y', 'FontSize', 18);
ylabel('Probability Density', 'FontSize', 18);
title('PDF of Normalized w_y Values (One Flow Field)', 'FontSize', 20);
set(gca, 'FontSize', 16);

% --- 3. Reference lines at x = 1.5 and x = 1.2 ----------------------------

% Vertical dashed lines
yl = ylim;
plot([1.5 1.5], yl, 'r--', 'LineWidth', 1.8);
plot([1.2 1.2], yl, 'r--', 'LineWidth', 1.8);

text(1.5, yl(2), '  x = 1.5', ...
     'Color','r', 'FontSize', 16, ...
     'VerticalAlignment','top', 'HorizontalAlignment','left');
text(1.2, yl(2), '  x = 1.2', ...
     'Color','r', 'FontSize', 16, ...
     'VerticalAlignment','top', 'HorizontalAlignment','right');

% --- 4. Shade Area Between x = 1.2 and x = 1.5 ----------------------------

% Find bin centers and counts between x = 1.2 and 1.5
inRange = centers >= 1.2 & centers <= 1.5;
x_fill = [centers(inRange), fliplr(centers(inRange))];
y_fill = [counts(inRange), zeros(1, sum(inRange))];

% Fill the area
fill(x_fill, y_fill, [0.3 0.7 1], 'FaceAlpha', 0.4, 'EdgeColor', 'none');

% --- 5. Estimate and Display Area (Probability) ---------------------------

dx = diff(edges);             % Bin widths
areaProb = sum(counts(inRange) .* dx(inRange));  % Approx. area

% Display probability as text
text(1.35, max(counts)*0.8, ...
     sprintf('P(1.2 ≤ x ≤ 1.5) ≈ %.4f', areaProb), ...
     'FontSize', 16, 'Color', 'b');





