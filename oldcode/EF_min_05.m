function E_min = EF_min_05(G)
    % Constants
    y = sqrt(0.5);  % ≈ 0.7071
    c = sqrt(0.5);  % ≈ 0.7071

    % Compute u that minimizes E
    u = (-0.7071 + sqrt(8*G.^2 - 11.3137*G + 8.5)) / 2;

    % Compute E using the original formula
    numerator = G ./ u;
    term = 1 + ((y^2 + u.^2 - 2*y*u + c^2) ./ G.^2) - (2*c ./ G);
    E_min = numerator .* term.^(3/4);
end
