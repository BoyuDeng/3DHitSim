function E_min = EF_min_01(G)
    % Constants
    y = sqrt(0.1);  % ≈ 0.3162
    c = sqrt(0.9);  % ≈ 0.9487

    % Compute u that minimizes E
    % Derived from: u^2 + y*u = 2G^2 - 4cG + 2(y^2 + c^2)
    % => u = (-y + sqrt(8G^2 - 16cG + 8(y^2 + c^2) + y^2)) / 2
    u = (-0.3162 + sqrt(8*G.^2 - 15.1792*G + 8.8016)) / 2;

    % Compute E using the original energy function
    numerator = G ./ u;
    term = 1 + ((0.1 + u.^2 - 2*0.3162*u + 0.9) ./ G.^2) - (2*0.9487 ./ G);
    E_min = numerator .* term.^(3/4);
end
