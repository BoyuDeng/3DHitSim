
dot_avgstotal=[];
updraft_avgstotal=[];
crosswind_avgstotal=[];
ywind_avgstotal=[];
mag_avgstotal = [];
Umag_avgstotal = [];
Wvalsall = [];
% Preallocate containers to collect all samples across all n and i
all_updraft = [];
all_ywind = [];
all_crosswind = [];

%%

numFields = length(results);        % e.g., 5
numSamples = length(results{1});    % e.g., 13 or so

dot_avgs = zeros(numFields, numSamples);  % Store dot_avg for each sample in each field
updraft_avgs = zeros(numFields, numSamples);
Vrms = zeros(numFields, numSamples);
 crosswind_avgs = zeros(numFields, numSamples);
 ywind_avgs = zeros(numFields, numSamples);
 mag_avgs = zeros(numFields, numSamples);
  Umag_avgs = zeros(numFields, numSamples);
Wvals = zeros(numFields, numSamples);  % 5 fields × 13 values


for n = 1:numFields

    fn = mod(n, 15);
    if fn == 0
       fn = 15;
    end



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
        coeffs = tempresults.coeffs;
        W = tempresults.W;

        Wvals(n, i) = W/U;

        X = X14unlim(t, coeffs, W);

        Velocity_fields = get_vel(uField, vField, wField, 64, X) ./ U;

        Vs = V14unlim(t, coeffs, W) ./ U;
        norm = sqrt(sum(Vs.^2,1));
        Vunit = Vs./ norm;

        Wnorm = sqrt(sum(Velocity_fields.^2,1));
        Wunit = Velocity_fields./ Wnorm;


        Wmag_time  = Wnorm;
        Umag_time = norm;
              
         



        dot_all = sum(Wunit .* Vunit,1);  % 1×T
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
        Wmag_mean  = mean(Wmag_time);
        Umag_mean = mean(Umag_time);
        
        dot_avgs(n, i) = dot_avg;
        updraft_avgs(n,i) = updraft_avg;
        crosswind_avgs(n,i) = crosswind_avg;
        ywind_avgs(n,i)     = ywind_avg;
        mag_avgs(n,i) = Wmag_mean;
        Umag_avgs(n,i) = Umag_mean;


        all_updraft   = [all_updraft, updraft_all];
        all_ywind     = [all_ywind, ywind_all];
        all_crosswind = [all_crosswind, crosswind_all];
    end
end

dot_avgstotal = [dot_avgstotal; dot_avgs];
updraft_avgstotal = [updraft_avgstotal; updraft_avgs];
crosswind_avgstotal = [crosswind_avgstotal; crosswind_avgs];
ywind_avgstotal   = [ywind_avgstotal; ywind_avgs]; 
mag_avgstotal = [mag_avgstotal; mag_avgs];
Umag_avgstotal = [Umag_avgstotal; Umag_avgs];


Wvalsall = [Wvalsall; Wvals];
% Collect all values




%%
% Means
mean_dot_per_field = mean(dot_avgstotal, 1);
mean_updraft = mean(updraft_avgstotal, 1);
mean_crosswind = mean(crosswind_avgstotal, 1);
mean_ywind = mean(ywind_avgstotal,  1);
meanWall = mean(Wvals, 1);
mean_mag  = mean(mag_avgstotal ,1);
mean_Umag  = mean(Umag_avgstotal ,1);

% Standard Errors
std_dot_per_field = std(dot_avgstotal, 0, 1) / sqrt(size(dot_avgstotal, 1));
std_updraft = std(updraft_avgstotal, 0, 1) / sqrt(size(updraft_avgstotal, 1));
std_crosswind = std(crosswind_avgstotal, 0, 1) / sqrt(size(crosswind_avgstotal, 1));
std_ywind = std(ywind_avgstotal, 0, 1) / sqrt(size(ywind_avgstotal,1));
std_wval = std(Wvalsall, 0, 1) / sqrt(size(Wvalsall,1));
std_mag   = std(mag_avgstotal ,0,1) / sqrt(size(mag_avgstotal ,1));
std_Umag   = std(Umag_avgstotal ,0,1) / sqrt(size(Umag_avgstotal ,1));
% Common font sizes
axisFont = 16;
labelFont = 18;
titleFont = 18;
legendFont = 16;

% Plot 1: <u · w> (Dot Product)
figure;
errorbar(Gall, mean_dot_per_field, std_dot_per_field, 'o-', 'LineWidth', 1.5);
set(gca, 'XScale', 'log', 'FontSize', axisFont);
xlabel('G', 'FontSize', labelFont);
ylabel('<W_{unit}· U_{unit}>', 'FontSize', labelFont);
title('Alignment of Flow and Trajectory', 'FontSize', titleFont);
grid on;

% Plot 2: Updraft
figure;
h1 = errorbar(Gall, mean_updraft, std_updraft, 's-', 'LineWidth', 1.5);
h2 = yline(-0.07, '--r', 'LineWidth', 1.5);
set(gca, 'XScale', 'log', 'FontSize', axisFont);
xlabel('G', 'FontSize', labelFont);
ylabel('<W_{z}>', 'FontSize', labelFont);
title('Mean Updraft', 'FontSize', titleFont);
legend([h1 h2], {'Mean Updraft', 'Mean Updraft for SLT'}, ...
       'Location', 'best', 'FontSize', legendFont);
grid on;

% Plot 3: Crosswind
figure;
h1 = errorbar(Gall, mean_crosswind, std_crosswind, 'd-', 'LineWidth', 1.5); hold on;
h2 = yline(0.78, '--r', 'LineWidth', 1.5);
set(gca, 'XScale', 'log', 'FontSize', axisFont);
xlabel('G', 'FontSize', labelFont);
ylabel('<|W_{x}|>', 'FontSize', labelFont);
title('Mean Crosswind', 'FontSize', titleFont);
legend([h1 h2], {'Mean Crosswind(CW)', 'Mean CW for Straight Line Trajectory'}, ...
       'Location', 'best', 'FontSize', legendFont);
grid on;

% Plot 4: ywind (Tailwind)
figure;
errorbar(Gall, mean_ywind, std_ywind, '^-', 'LineWidth', 1.5);
set(gca, 'XScale', 'log', 'FontSize', axisFont);
xlabel('G', 'FontSize', labelFont);
ylabel('<W_{y}>', 'FontSize', labelFont);
title('Mean Tailwind', 'FontSize', titleFont);
grid on;

% Plot 5: W value (U_0)
figure;
errorbar(Gall, meanWall, std_wval, '^-', 'LineWidth', 1.5); hold on;
% Add y = sqrt(2)*x line
xline = logspace(log10(min(Gall)), log10(max(Gall)), 100); % generate 100 log-spaced points
ylineqf = sqrt(2) * xline;
plot(xline, ylineqf, 'k-.', 'LineWidth', 1.5);


set(gca, 'XScale', 'log', 'FontSize', axisFont);
xlabel('G', 'FontSize', labelFont);
ylabel('<U_0>', 'FontSize', labelFont);
title('', 'FontSize', titleFont);
legend('<U_0>', 'U for QF model', 'Location', 'best', 'FontSize', legendFont);
grid on;


% Plot 6: Velocity-magnitude W and U
figure;
errorbar(Gall, mean_mag, std_mag, 'v-', 'LineWidth', 1.5); hold on;
errorbar(Gall, mean_Umag, std_Umag, 'x--', 'LineWidth', 1.5);

% Add y = sqrt(2)*x line
xline = logspace(log10(min(Gall)), log10(max(Gall)), 100); % generate 100 log-spaced points
ylineqf = sqrt(2) * xline;
plot(xline, ylineqf, 'k-.', 'LineWidth', 1.5);

set(gca, 'XScale', 'log', 'FontSize', axisFont);
xlabel('G', 'FontSize', labelFont);
ylabel('Velocity Magnitude', 'FontSize', labelFont);
title('Mean Velocity Magnitude', 'FontSize', titleFont);
legend('<W> (flow field)', '<U> (point-mass)', 'U for QF model', 'Location', 'best', 'FontSize', legendFont);
grid on;

%% PDF of Updraft  (w_z)
figure;
datasets = {all_updraft_lo, all_updraft_mi, all_updraft_hi};
labels = {'G 0.5 to 1', '1 to 6', '6 to 25'};
colors = [0 0.50 0.80; 0.00 0.60 0.00; 0.80 0.00 0.00];

for i = 1:3
    [counts, edges] = histcounts(datasets{i}, 'BinWidth', 0.05, 'Normalization', 'pdf');
    centers = edges(1:end-1) + diff(edges)/2;
    stairs(centers, counts, 'Color', colors(i,:), 'LineWidth', 1.8); hold on;
end
xlabel('w_z (Updraft)', 'FontSize', 18);
ylabel('Probability Density', 'FontSize', 18);
title('PDF of w_z (Updraft)', 'FontSize', 20);
legend(labels, 'FontSize', 14);
grid on;
set(gca, 'FontSize', 16);

%% PDF of Tailwind  (w_y)
figure;
datasets = {all_ywind_lo, all_ywind_mi, all_ywind_hi};

for i = 1:3
    [counts, edges] = histcounts(datasets{i}, 'BinWidth', 0.05, 'Normalization', 'pdf');
    centers = edges(1:end-1) + diff(edges)/2;
    stairs(centers, counts, 'Color', colors(i,:), 'LineWidth', 1.8); hold on;
end
xlabel('w_y (Tailwind)', 'FontSize', 18);
ylabel('Probability Density', 'FontSize', 18);
title('PDF of w_y (Tailwind)', 'FontSize', 20);
legend(labels, 'FontSize', 14);
grid on;
set(gca, 'FontSize', 16);


%% PDF of Crosswind  (|w_x|)
figure;
datasets = {all_crosswind_lo, all_crosswind_mi, all_crosswind_hi};

for i = 1:3
    [counts, edges] = histcounts(datasets{i}, 'BinWidth', 0.05, 'Normalization', 'pdf');
    centers = edges(1:end-1) + diff(edges)/2;
    stairs(centers, counts, 'Color', colors(i,:), 'LineWidth', 1.8); hold on;
end
xlabel('|w_x| (Crosswind)', 'FontSize', 18);
ylabel('Probability Density', 'FontSize', 18);
title('PDF of |w_x| (Crosswind)', 'FontSize', 20);
legend(labels, 'FontSize', 14);
grid on;
set(gca, 'FontSize', 16);









