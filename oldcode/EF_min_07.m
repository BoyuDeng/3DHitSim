function E_min = EF_min_07(G)
    % Constants
    y = sqrt(0.7);  % ≈ 0.8367
    c = sqrt(0.3);  % ≈ 0.5477

    % Compute u that minimizes E
    u = (-0.8367 + sqrt(8*G.^2 - 8.7632*G + 8.6997)) / 2;

    % Compute E using the given formula
    numerator = G ./ u;
    term = 1 + ((y^2 + u.^2 - 2*y*u + c^2) ./ G.^2) - (2*c ./ G);
    E_min = numerator .* term.^(3/4);
end
