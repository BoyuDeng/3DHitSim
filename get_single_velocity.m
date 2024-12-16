function vel = get_single_velocity(pos, U, V, W, sizeo)
    % Initialize the velocity output array
    % Ensure position stays in [0, 1] box with periodic boundary conditions
    pos = mod(pos, 1); 
    vel = zeros(3,1);
    
    x = linspace(0,1,sizeo);
    y = x;
    z = x;
    xm =linspace(x(2)/2,x(end)+x(2)/2,sizeo);
    ym = xm;
    zm = xm;

    vect = get_particle_cell_Location(pos, sizeo);
    i0 = vect(1);
    j0 = vect(2);
    k0 = vect(3);
    imaxo = sizeo;
    jmaxo = sizeo;
    kmaxo = sizeo;
    imino = 1;
    jmino = 1;
    kmino = 1;


    % Interpolate U velocity ------------------------------
    % Find right i index for U
    i = max(min(imaxo - 1, i0), imino);
    while (pos(1) - x(i) < 0.0 && i > imino)
        i = i - 1;
    end
    while (pos(1) - x(i + 1) >= 0.0 && i + 1 < imaxo)
        i = i + 1;
    end

    % Find right j index for U
    j = max(min(jmaxo - 1, j0), jmino);
    while (pos(2) - ym(j) < 0.0 && j > jmino)
        j = j - 1;
    end
    while (pos(2) - ym(j + 1) >= 0.0 && j + 1 < jmaxo)
        j = j + 1;
    end

    % Find right k index for U
    k = max(min(kmaxo - 1, k0), kmino);
    while (pos(3) - zm(k) < 0.0 && k > kmino)
        k = k - 1;
    end
    while (pos(3) - zm(k + 1) >= 0.0 && k + 1 < kmaxo)
        k = k + 1;
    end

    % Prepare tri-linear interpolation coefficients for U
    wx1 = (pos(1) - x(i)) / (x(i + 1) - x(i)); wx2 = 1.0 - wx1;
    wy1 = (pos(2) - ym(j)) / (ym(j + 1) - ym(j)); wy2 = 1.0 - wy1;
    wz1 = (pos(3) - zm(k)) / (zm(k + 1) - zm(k)); wz2 = 1.0 - wz1;

    % Tri-linear interpolation of U
    vel(1) = wz1 * (wy1 * (wx1 * U(i + 1, j + 1, k + 1) + wx2 * U(i, j + 1, k + 1)) + ...
                    wy2 * (wx1 * U(i + 1, j, k + 1) + wx2 * U(i, j, k + 1))) + ...
             wz2 * (wy1 * (wx1 * U(i + 1, j + 1, k) + wx2 * U(i, j + 1, k)) + ...
                    wy2 * (wx1 * U(i + 1, j, k) + wx2 * U(i, j, k)));

    % Interpolate V velocity ------------------------------
    % Find right i index for V
    i = max(min(imaxo - 1, i0), imino);
    while (pos(1) - xm(i) < 0.0 && i > imino)
        i = i - 1;
    end
    while (pos(1) - xm(i + 1) >= 0.0 && i + 1 < imaxo)
        i = i + 1;
    end

    % Find right j index for V
    j = max(min(jmaxo - 1, j0), jmino);
    while (pos(2) - y(j) < 0.0 && j > jmino)
        j = j - 1;
    end
    while (pos(2) - y(j + 1) >= 0.0 && j + 1 < jmaxo)
        j = j + 1;
    end

    % Find right k index for V
    k = max(min(kmaxo - 1, k0), kmino);
    while (pos(3) - zm(k) < 0.0 && k > kmino)
        k = k - 1;
    end
    while (pos(3) - zm(k + 1) >= 0.0 && k + 1 < kmaxo)
        k = k + 1;
    end

    % Prepare tri-linear interpolation coefficients for V
    wx1 = (pos(1) - xm(i)) / (xm(i + 1) - xm(i)); wx2 = 1.0 - wx1;
    wy1 = (pos(2) - y(j)) / (y(j + 1) - y(j)); wy2 = 1.0 - wy1;
    wz1 = (pos(3) - zm(k)) / (zm(k + 1) - zm(k)); wz2 = 1.0 - wz1;

    % Tri-linear interpolation of V
    vel(2) = wz1 * (wy1 * (wx1 * V(i + 1, j + 1, k + 1) + wx2 * V(i, j + 1, k + 1)) + ...
                    wy2 * (wx1 * V(i + 1, j, k + 1) + wx2 * V(i, j, k + 1))) + ...
             wz2 * (wy1 * (wx1 * V(i + 1, j + 1, k) + wx2 * V(i, j + 1, k)) + ...
                    wy2 * (wx1 * V(i + 1, j, k) + wx2 * V(i, j, k)));

    % Interpolate W velocity ------------------------------
    % Find right i index for W
    i = max(min(imaxo - 1, i0), imino);
    while (pos(1) - xm(i) < 0.0 && i > imino)
        i = i - 1;
    end
    while (pos(1) - xm(i + 1) >= 0.0 && i + 1 < imaxo)
        i = i + 1;
    end

    % Find right j index for W
    j = max(min(jmaxo - 1, j0), jmino);
    while (pos(2) - ym(j) < 0.0 && j > jmino)
        j = j - 1;
    end
    while (pos(2) - ym(j + 1) >= 0.0 && j + 1 < jmaxo)
        j = j + 1;
    end

    % Find right k index for W
    k = max(min(kmaxo - 1, k0), kmino);
    while (pos(3) - z(k) < 0.0 && k > kmino)
        k = k - 1;
    end
    while (pos(3) - z(k + 1) >= 0.0 && k + 1 < kmaxo)
        k = k + 1;
    end

    % Prepare tri-linear interpolation coefficients for W
    wx1 = (pos(1) - xm(i)) / (xm(i + 1) - xm(i)); wx2 = 1.0 - wx1;
    wy1 = (pos(2) - ym(j)) / (ym(j + 1) - ym(j)); wy2 = 1.0 - wy1;
    wz1 = (pos(3) - z(k)) / (z(k + 1) - z(k)); wz2 = 1.0 - wz1;

    % Tri-linear interpolation of W
    vel(3) = wz1 * (wy1 * (wx1 * W(i + 1, j + 1, k + 1) + wx2 * W(i, j + 1, k + 1)) + ...
                    wy2 * (wx1 * W(i + 1, j, k + 1) + wx2 * W(i, j, k + 1))) + ...
             wz2 * (wy1 * (wx1 * W(i + 1, j + 1, k) + wx2 * W(i, j + 1, k)) + ...
                    wy2 * (wx1 * W(i + 1, j, k) + wx2 * W(i, j, k)));
end
