function [E_min, u] = EF_min_general(y, c,s, G)
    % EF_min_general: Computes minimum E for given y, c, and G
    %
    % Inputs:
    %   y - scalar (sqrt of y^2)
    %   c - scalar (sqrt of c^2)
    %   G - scalar or vector of G values
    %
    % Output:
    %   E_min - minimum energy values corresponding to each G

    % Compute constants
    y2 = y^2;
    c2 = c^2;
    s2 = s^2;

    % Solve for u_crit
    u = (-y + sqrt(8*G.^2 - 16*c*G + 8*(y^2 + c^2+ s^2) + y^2)) / 2;


    % Energy function
    term = 1 + ((y2 + u.^2 - 2*y.*u + c2 + s2) ./ G.^2) - (2*c ./ G);
    E_min = (G ./ u) .* term.^(3/4);
end
