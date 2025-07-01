function E_min = EF_min_y25_c04(G)
    % Constants
    y = 2.5;
    c = 0.4;

    % Solve for u_crit
    u = (-2.5 + sqrt(8*G.^2 - 6.4*G + 25.64)) / 2;

    % Compute E using original energy function
    numerator = G ./ u;
    term = 1 + ((y^2 + u.^2 - 2*y*u + c^2) ./ G.^2) - (2*c ./ G);
    E_min = numerator .* term.^(3/4);
end
