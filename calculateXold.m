function x = calculateXold(t, a, b, c, tf, W, a0, c0, b0)
    % Determine the number of terms in the sums
    n = length(a);

    % Initialize the components of the vector
    x1 = a0;  % Assuming a0 is the first element of vector a
    x2 = b0 + W * t ;
    x3 = c0;          % Assuming c0 is the first element of vector c

    % Compute the sums for each component
    for i = 1:n
        x1 = x1 + a(i) * sin(i * pi * t / tf);
        x2 = x2 + b(i) * sin(i * pi * t / tf);
        x3 = x3 + c(i) * sin(i * pi * t / tf);
    end

    % Combine the components into a column vector
    x = [x1; x2; x3];
end