


function X = calculateX(t, coeffs, W, StartLoc)
    % Determine the number of terms in the sums
    n = length(coeffs) / 3;
    
    % Extract coefficients for a, b, and c from the single coeffs vector
    a = coeffs(1:n);
    b = coeffs(n+1:2*n);
    c = coeffs(2*n+1:end);
    
    % Set odd-indexed elements of b to zero
    b(1:2:end) = 0;
    
    % Initialize the components of the vector
    a0 = StartLoc(1);
    b0 = 0.1;
    c0 = StartLoc(2);  % Assuming c0 is initially 0
    
    % Initialize the result matrix X
    X = zeros(3, length(t));
    
    % Compute the sums for each component
    for i = 1:n
        X(1, :) = X(1, :) + a(i) * sin(i * pi * t / (2*t(end)));
        X(2, :) = X(2, :) + b(i) * sin(i * pi * t / t(end));
        X(3, :) = X(3, :) + c(i) * sin(i * pi * t / (2*t(end)));
    end
    
    % Add the initial positions and linear term for x2
    X(1, :) = X(1, :) + a0;
    X(2, :) = X(2, :) + b0 + W * t;
    X(3, :) = X(3, :) + c0;
end

