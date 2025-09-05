function E_min = EF_min_02(G)
    % Constants
    y = sqrt(0.2);  % ≈ 0.4472
    c = sqrt(0.8);  % ≈ 0.8944

    % Compute u that minimizes E
    u = (-0.4472 + sqrt(8*G.^2 - 14.3104*G + 8.9999)) / 2;

    % Compute E using the original formula
    numerator = G ./ u;
    term = 1 + ((y^2 + u.^2 - 2*y*u + c^2) ./ G.^2) - (2*c ./ G);
    E_min = numerator .* term.^(3/4);
end
