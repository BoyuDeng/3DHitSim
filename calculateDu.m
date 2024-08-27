function V = calculateDu(t, coeffs)
    % Determine the number of terms in the sums
    n = length(coeffs) / 3;
    
    % Extract coefficients for a, b, and c from the single coeffs vector
    a = coeffs(1:n);
    b = coeffs(n+1:2*n);
    c = coeffs(2*n+1:end);
    
    % Set odd-indexed elements of b to zero
    b(1:2:end) = 0;
    
    % Initialize the result matrix V
    V = zeros(3, length(t));
    
    % Compute the derivative for each component
    for i = 1:n
        V(1, :) = V(1, :) + a(i) * (i * pi / (2*t(end)))^2 * -sin(i * pi * t / (2*t(end)));
        V(2, :) = V(2, :) + b(i) * (i * pi / t(end))^2 * -sin(i * pi * t / t(end));
        V(3, :) = V(3, :) + c(i) * (i * pi / (2*t(end)))^2 * -sin(i * pi * t / (2*t(end)));
    end
    
end