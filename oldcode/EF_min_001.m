function E_min = EF_min_001(G)
    % Constants
    y = 0.1;
    c = sqrt(0.99);  % ≈ 0.9949

    % Compute u_crit
    u = (-y + sqrt(8*G.^2 - 4*4*c*G + 8*(y^2 + c^2))) / 2;

    % Compute E_min
    numerator = G ./ u;
    term = 1 + ((y^2 + u.^2 - 2*y*u + c^2) ./ G.^2) - (2*c ./ G);
    E_min = numerator .* term.^(3/4);
end
