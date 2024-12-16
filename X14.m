function X = X14(t, coeffs, W)
    % Validate inputs
    if nargin < 3
        error('X14 requires three inputs: t, coeffs, and W.');
    end
    if length(coeffs) < 3
        error('Coeffs must have at least 3 elements.');
    end
    if mod(length(coeffs) - 2, 3) ~= 0
        error('Coeffs length must satisfy length(coeffs) = 3n + 2.');
    end
    if isempty(t) || ~isnumeric(t)
        error('Time vector t must be a non-empty numeric array.');
    end
    if ~isscalar(W) || ~isnumeric(W)
        error('W must be a scalar numeric value.');
    end

    % Determine the number of terms in the sums
    n = (length(coeffs) - 2) / 3;

    % Extract coefficients for a, b, and c from the single coeffs vector
    a = coeffs(1:n);
    b = coeffs(n+1:2*n);
    c = coeffs(2*n+1:end-2);

    % Ensure coefficients are non-empty and valid
    if isempty(a) || isempty(b) || isempty(c)
        error('Invalid coefficients: a, b, or c is empty.');
    end

    % Set odd-indexed elements of b to zero
    %b(1:2:end) = 0;
    b(:) = 0;

    % Initialize constants
    last = coeffs(end-1);
    a0 = last + 0.3;
    b0 = 0.1;
    c0 = coeffs(end) + 0.7;

    % Initialize the result matrix X
    X = zeros(3, length(t));

    % Compute the sums for each component
    for i = 1:n
        X(1, :) = X(1, :) + a(i) * sin(i * pi * t / (2 * t(end)));
        X(2, :) = X(2, :) + b(i) * sin(i * pi * t / t(end));
        X(3, :) = X(3, :) + c(i) * sin(i * pi * t / (2 * t(end)));
    end

    % Add the initial positions and linear term for x2
    X(1, :) = X(1, :) + a0;
    X(2, :) = X(2, :) + b0 + W * t;
    X(3, :) = X(3, :) + c0;

    % Validate that X is within [0, 1]
    if any(X(:) < 0) || any(X(:) > 1)
        warning('X contains values outside the range [0, 1].');
    end
end
