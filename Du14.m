function V = Du14(t, coeffs)
    % Determine the number of terms in the sums
    n = (length(coeffs) - 2)*(2 / 5);
    
    % Extract coefficients for a, b, and c from the single coeffs vector
    a = coeffs(1:n);
    b = zeros(length(a),1);
    b(2:2:end) = coeffs(n+1:3*n/2);
    c = coeffs((3*n/2)+1:end-2);
    
    
    
    % Initialize the result matrix V
    V = zeros(3, length(t));
    
    % Compute the derivative for each component
    for i = 1:n
        V(1, :) = V(1, :) + a(i) * (i * pi / (2*t(end)))^2 * -sin(i * pi * t / (2*t(end)));
        V(2, :) = V(2, :) + b(i) * (i * pi / t(end))^2 * -sin(i * pi * t / t(end));
        V(3, :) = V(3, :) + c(i) * (i * pi / (2*t(end)))^2 * -sin(i * pi * t / (2*t(end)));
    end
    
end