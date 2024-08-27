function vectors = get_vectors_at_t(uField, vField, wField, t_index)
    % Validate the time index
    if t_index < 1 || t_index > numel(uField)
        error('Invalid time index.');
    end
    
    % Extract the velocity fields at the given time step
    u = uField{t_index};
    u = u{1};% uField is a cell array, access its contents with {}
    v = vField{t_index};
    v = v{1};
    w = wField{t_index};
    w = w{1};
    
    % Initialize the array to hold vectors
    vectors = zeros(64, 64, 64, 3);
    
    % Iterate through each point in the mesh
    for i = 1:64
        for j = 1:64
            for k = 1:64
                % Extract the velocity vector at this point
                vector = [u(i, j, k), v(i, j, k), w(i, j, k)];
                
                % Store the vector
                vectors(i, j, k, :) = vector;
            end
        end
    end
end