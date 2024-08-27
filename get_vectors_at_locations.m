function vectors = get_vectors_at_locations(uField, vField, wField, cellIndices)
    % Validate the time steps
    numTimeSteps = size(cellIndices, 2);
    if numTimeSteps > numel(uField)
        error('Number of time steps exceeds the available data.');
    end
    
    % Initialize the array to hold velocity vectors
    vectors = zeros(3, numTimeSteps);
    
    % Iterate through each time step
    for t = 1:numTimeSteps
        % Validate the spatial indices for the current time step
        i = cellIndices(1, t);
        j = cellIndices(2, t);
        k = cellIndices(3, t);
        
        if i < 1 || i > 64 || j < 1 || j > 64 || k < 1 || k > 64
            error('Invalid spatial indices at time step %d. i, j, and k must be within the range 1 to 64.', t);
        end
        
        % Extract the velocity fields at the current time step
        u = uField{t};
        u = u{1}; % uField is a cell array, access its contents with {}
        v = vField{t};
        v = v{1};
        w = wField{t};
        w = w{1};
        
        % Extract the velocity vector at the specified location
        vectors(:, t) = [u(i, j, k), v(i, j, k), w(i, j, k)];
    end
end