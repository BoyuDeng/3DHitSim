% --- user parameters -----------------------------------------------------
modes      = [10 20 ];                       % retained-mode counts
resultSets = {resultsU8_10, resultsU8_20, ...
              };        % 4 optimisation levels
nFields    = 5;                                   % outer cells to average
realIdx    =12;                                  % realisation index
components = {'u','v','w'};                       % which velocity fields
% -------------------------------------------------------------------------

nComp = numel(components);
figure;

% Font settings
axisFont = 14;
labelFont = 16;
titleFont = 18;
legendFont = 14;

for c = 1:nComp
    subplot(1, 3, c); hold on
    for k = 1:numel(modes)
        m   = modes(k);
        Ec  = zeros(nFields, m);  % energy per field, per mode index

        for f = 1:nFields
            coeffs = resultSets{k}{f}{realIdx}.coeffs;
            offset = (c-1)*m;  % 0 for u, m for v, 2m for w
            Ec(f,:) = coeffs(1+offset : m+offset).^2;
        end

        Ec_mean = mean(Ec, 1);  % average over the 5 fields
        h(k) = plot(1:m, Ec_mean, '-o', ...
                    'DisplayName', sprintf('%d', m), ...
                    'LineWidth', 1.5);
    end

    ax = gca;
    set(ax, 'YScale', 'log', 'Box', 'on', ...
        'GridLineStyle', ':', 'FontSize', axisFont);

    xlabel('Mode number', 'FontSize', labelFont);
    ylabel(sprintf('%s-energy', components{c}), 'FontSize', labelFont);

    if c == 1
        legend('Location', 'southwest', 'FontSize', legendFont);
    end

    grid on;
end

% Add overall title
sgtitle('Spectrum (mean of 5 trials, G = 2)', 'FontSize', titleFont);
% Get current figure handle
f = gcf;

% Set figure size: [left, bottom, width, height]
% Increase width (e.g., 1200 px) to stretch horizontally
% --- Save moderately wide PDF version ---
f.Position = [100, 100, 900, 400];   % Slightly wider than compact
exportgraphics(f, 'spectrum_plot_G20.pdf', 'ContentType', 'vector');
