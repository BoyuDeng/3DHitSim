function E_min = EF_min_08(G)
    % Constants
    y = sqrt(0.8);  % ≈ 0.8944
    c = sqrt(0.2);  % ≈ 0.4472

    % Compute u that minimizes E
    u = (-0.8944 + sqrt(8*G.^2 - 7.1552*G + 8.7994)) / 2;

    % Compute E using the given formula
    numerator = G ./ u;
    term = 1 + ((y^2 + u.^2 - 2*y*u + c^2) ./ G.^2) - (2*c ./ G);
    E_min = numerator .* term.^(3/4);
end
