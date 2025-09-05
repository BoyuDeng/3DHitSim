function L = calculate_scale(uField, dx)
    % Function to calculate the integral length scale from velocity data
    %
    % Inputs:
    %   u  - Vector of velocity data
    %   dx - Spatial (or temporal) increment between samples
    %
    % Output:
    %   L  - Integral length scale
    u = zeros(650,1);
    % Ensure input is a column vector
    for t  = 1:650
        a = uField{t};
        a = a{1};
        u(t) = a(1);
    end
    
    % Remove the mean from the velocity signal
    u_mean = mean(u);
    u_fluct = u - u_mean;
    
    % Calculate the autocorrelation of the velocity field
    [autocorr_u, lags] = xcorr(u_fluct, 'biased');
    
    % Normalize the autocorrelation
    autocorr_u = autocorr_u / max(autocorr_u);
    
    % Find the positive lags
    positive_lags = lags(lags >= 0);
    autocorr_u_positive = autocorr_u(lags >= 0);
    
    % Estimate the integral length scale (L)
    L = trapz(positive_lags * dx, autocorr_u_positive);
    
    % Display the result
    disp(['Integral Length Scale: ', num2str(L)]);
end
