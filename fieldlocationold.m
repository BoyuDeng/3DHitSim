function vector = get_vector_at_location(uField, vField, wField, t_index, i, j, k)
    % Validate the time index
    if t_index < 1 || t_index > numel(uField)
        error('Invalid time index.');
    end

    % Validate the spatial indices
    if i < 1 || i > 64 || j < 1 || j > 64 || k < 1 || k > 64
        error('Invalid spatial indices. i, j, and k must be within the range 1 to 64.');
    end
    
    % Extract the velocity fields at the given time step
    u = uField{t_index};
    u = u{1};% uField is a cell array, access its contents with {}
    v = vField{t_index};
    v = v{1};
    w = wField{t_index};
    w = w{1};
    
    % Extract the velocity vector at the specified location
    vector = [u(i, j, k), v(i, j, k), w(i, j, k)];
end


