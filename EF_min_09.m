function E_min = EF_min_09(G)
    % Constants
    y = sqrt(0.9);
    c = sqrt(0.1);
    % Compute u that minimizes E
    u = (-0.9486 + sqrt(8*G.^2 - 5.0592*G + 8.8999)) / 2;

    % Compute E using the given formula
    numerator = G ./ u;
    term = 1 + ((y^2 + u.^2 - 2*y*u + c^2) ./ G.^2) - (2*c ./ G);
    E_min = numerator .* term.^(3/4);
end
