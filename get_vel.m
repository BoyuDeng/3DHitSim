function vectors = get_vel(uField, vField, wField, sizeo, X)
    % Validate the time steps
    numTimeSteps = size(X, 2);
    
    if numTimeSteps > numel(uField)
        error('Number of time steps exceeds the available data.');
    end
    
    % Initialize the array to hold velocity vectors
    vectors = zeros(3, numTimeSteps);
    
    % Iterate through each time step
    for t = 1:numTimeSteps
        u = uField{t};
        u = u{1}; 
        v = vField{t};
        v = v{1};
        w = wField{t};
        w = w{1};
        vectors(:, t) = get_single_velocity(X(:,t),u,v,w,sizeo);
    end
end