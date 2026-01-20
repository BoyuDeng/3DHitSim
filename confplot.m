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

Gall = [ 1, 25 , 5, 8 , 11, 15];

Gall = [1.5, 2 , 3, 4];

Gall = [5,6,8,10];

Gall = [1.5, 2 , 3, 4, 5,6,8,10, 12, 16, 20 25];

% [Gall_sorted, idx] = sort(Gallunsort);
% results = resultsunsort(idx);

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



%  % Define the file path dynamically
% file_path = sprintf('/Users/boyi/Simulations/hit/output_fields%dU=8', 100);
% 
% % Define indices to read (every 3rd file: 1, 4, 7, ...)
% indices = start_num:1:end_num;
% numFiles = numel(indices);
% variables_list = cell(numFiles, 1);
% 
% % Read selected field files
% for k = 1:numFiles
%     i = indices(k);
%     filename = fullfile(file_path, sprintf('field_%05d', i));
%     [~, ~, ~, ~, variables] = read_field_file(filename);
%     variables_list{k} = variables;
% end
% 
% % Extract u, v, w fields
% uField = cell(numFiles, 1);
% vField = cell(numFiles, 1);
% wField = cell(numFiles, 1);
% 
% for k = 1:numFiles
%     data = variables_list{k};
%     uField{k} = data(1);
%     vField{k} = data(2);
%     wField{k} = data(3);
% end

% Define time vector to match sampled indices
t = (indices - start_num) * dt;
% Define the file path dynamically
file_path = sprintf('/Users/boyi/Simulations/hit/output_fields%dU=8', 100);

% Define indices to read (every 3rd file: 1, 4, 7, ...)
indices = start_num:1:end_num;
numFiles = numel(indices);
variables_list = cell(numFiles, 1);

% Read selected field files
for k = 1:numFiles
    i = indices(k);
    filename = fullfile(file_path, sprintf('field_%05d', i));
    [~, ~, ~, ~, variables] = read_field_file(filename);
    variables_list{k} = variables;
end

% Extract u, v, w fields
uField = cell(numFiles, 1);
vField = cell(numFiles, 1);
wField = cell(numFiles, 1);

for k = 1:numFiles
    data = variables_list{k};
    uField{k} = data(1);
    vField{k} = data(2);
    wField{k} = data(3);
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

% Use st_Gall for all st-series
x_st = st_Gall;

% Plot 1: <u · w> (Dot Product)
figure;

h_st  = errorbar(x_st, st_mean_dot_per_field, st_std_dot_per_field, ...
    'o', 'MarkerFaceColor', 'auto', 'LineWidth', 1.5, 'MarkerSize', 12);

set(gca, 'XScale', 'log', 'FontSize', axisFont);
xlabel('G', 'FontSize', labelFont);
ylabel('<W_{unit}\cdot U_{unit}>', 'FontSize', labelFont);
title('Alignment of Flow and Trajectory', 'FontSize', titleFont);
%legend([h_reg h_st], {'Regular','st'}, 'Location','best', 'FontSize', legendFont);

grid on;
ax = gca;
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
%axis([5^-2 50 0 1])



% Plot 2: Updraft
figure;
%h_reg = errorbar(Gall, mean_updraft, std_updraft, 's', 'LineWidth', 1.5); hold on;

h_st  = errorbar(x_st, st_mean_updraft, st_std_updraft, 's', 'MarkerFaceColor', 'auto', 'LineWidth', 1.5, 'MarkerSize', 12);
h_ref = yline(-0.07, '--r', 'LineWidth', 1.5);
set(gca, 'XScale', 'log', 'FontSize', axisFont);
xlabel('G', 'FontSize', labelFont);
ylabel('<W_{z}>', 'FontSize', labelFont);
title('Mean Updraft', 'FontSize', titleFont);
%legend([h_reg h_st h_ref], {'Regular','st','Mean Updraft for SLT'}, 'Location','best','FontSize',legendFont);
grid on;
ax = gca;
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
axis([5^-2 50 -0.2 1.2])


% Plot 3: Crosswind
figure;
%h_reg = errorbar(Gall, mean_crosswind, std_crosswind, 'd', 'LineWidth', 1.5); hold on;

h_st  = errorbar(x_st, st_mean_crosswind, st_std_crosswind, 'd', 'MarkerFaceColor', 'auto', 'LineWidth', 1.5, 'MarkerSize', 12);
h_ref = yline(0.78, '--r', 'LineWidth', 1.5);
set(gca, 'XScale', 'log', 'FontSize', axisFont);
xlabel('G', 'FontSize', labelFont);
ylabel('<|W_{x}|>', 'FontSize', labelFont);
title('Mean Crosswind', 'FontSize', titleFont);
grid on;
ax = gca;
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
axis([5^-2 50 0 0.9])

% Plot 4: ywind (Tailwind)
figure;
%h_reg = errorbar(Gall, mean_ywind, std_ywind, '^', 'LineWidth', 1.5); hold on;
h_st  = errorbar(x_st, st_mean_ywind, st_std_ywind, '^', 'MarkerFaceColor', 'auto', 'LineWidth', 1.5, 'MarkerSize', 12);
set(gca, 'XScale', 'log', 'FontSize', axisFont);
xlabel('G', 'FontSize', labelFont);
ylabel('<W_{y}>', 'FontSize', labelFont);
title('Mean Tailwind', 'FontSize', titleFont);
%legend([h_reg h_st], {'Regular','st'}, 'Location','best','FontSize',legendFont);
grid on;
ax = gca;
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
axis([5^-2 50 0 1.7])

% % Plot 5: W value (U_0)
% figure;
% h_reg = errorbar(Gall, meanWall, std_wval, '^', 'LineWidth', 1.5); hold on;
% h_st  = errorbar(x_st, st_meanWall, st_std_wval, '^', 'MarkerFaceColor', 'auto', 'LineWidth', 1.5, 'MarkerSize', 12);
% 
% % Add y = sqrt(2)*x line (based on regular Gall extent)
% x_qf = logspace(log10(min(Gall)), log10(max(Gall)), 100);
% y_qf = sqrt(2) * x_qf;
% h_qf = plot(x_qf, y_qf, 'k-.', 'LineWidth', 1.5);
% 
% set(gca, 'XScale', 'log', 'FontSize', axisFont);
% xlabel('G', 'FontSize', labelFont);
% ylabel('<U_0>', 'FontSize', labelFont);
% title('', 'FontSize', titleFont);
% legend([h_reg h_st h_qf], {'Regular','st','U for QF model'}, 'Location','best','FontSize',legendFont);
% grid on;

% Plot 6: Velocity-magnitude W and U
figure; hold on

h_W_st = errorbar(x_st, st_mean_mag,  st_std_mag, ...
    'v', 'MarkerFaceColor', 'auto', 'LineWidth', 1.5, 'MarkerSize', 12, ...
    'DisplayName','<W> (flow field)');

h_U_st = errorbar(x_st, st_mean_Umag, st_std_Umag, ...
    's', 'MarkerFaceColor', 'auto', 'LineWidth', 1.5, 'MarkerSize', 12, ...
    'DisplayName','<U> (point-mass)');

% Add y = sqrt(2)*x line (based on regular Gall extent)
x_qf = logspace(log10(min(Gall)), log10(max(Gall)), 100);
y_qf = sqrt(2) * x_qf;
h_qf = plot(x_qf, y_qf, 'k-.', 'LineWidth', 1.5, 'DisplayName','U for QF model');

set(gca, 'XScale', 'log', 'FontSize', axisFont);
xlabel('G', 'FontSize', labelFont);
ylabel('Velocity Magnitude', 'FontSize', labelFont);
title('Mean Velocity Magnitude', 'FontSize', titleFont);

legend('Location','best','FontSize',legendFont);  % uses DisplayName

grid on;
ax = gca;
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
axis([5^-2 50 -0.1 36])


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


%%

save('st_data.mat', ...
    'st_mean_dot_per_field','st_mean_updraft','st_mean_crosswind','st_mean_ywind','st_meanWall','st_mean_mag','st_mean_Umag', ...
    'st_std_dot_per_field','st_std_updraft','st_std_crosswind','st_std_ywind','st_std_wval','st_std_mag','st_std_Umag','st_Gall');

save('Sep24data.mat', 'resultsallU8_4_10x','resultsallU8_5_10x','resultsallU8_6_10x')


save('Dec3data.mat', 'resultsallU8_3x','resultsallU8_10x','resultsallU8_5x','resultsallU8_10x2')

