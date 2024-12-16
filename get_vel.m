function vectors = get_vel(uField, vField, wField, sizeo, X)
    % Validate inputs
    if nargin < 5
        error('get_vel requires five inputs: uField, vField, wField, sizeo, and X.');
    end
    if ~iscell(uField) || ~iscell(vField) || ~iscell(wField)
        error('uField, vField, and wField must be cell arrays.');
    end
    if size(X, 1) ~= 3
        error('X must have 3 rows corresponding to the x, y, z components.');
    end
    
    % Validate time steps
    numTimeSteps = size(X, 2);
    if numTimeSteps > numel(uField) || numTimeSteps > numel(vField) || numTimeSteps > numel(wField)
        error('Number of time steps exceeds the available velocity field data.');
    end
    
    % Initialize the array to hold velocity vectors
    vectors = zeros(3, numTimeSteps);
    
    % Iterate through each time step
    for t = 1:numTimeSteps
        % Extract velocity fields for the current time step
        try
            u = uField{t};
            v = vField{t};
            w = wField{t};
        catch
            error('Error accessing velocity fields at time step %d.', t);
        end

        % Validate extracted fields
        if ~iscell(u) || ~iscell(v) || ~iscell(w) || isempty(u{1}) || isempty(v{1}) || isempty(w{1})
            error('Velocity fields at time step %d are invalid or empty.', t);
        end
        
        % Extract the field arrays
        u = u{1};
        v = v{1};
        w = w{1};

        % Check field size matches `sizeo`
        if size(u, 1) ~= sizeo || size(u, 2) ~= sizeo || size(v, 1) ~= sizeo || size(v, 2) ~= sizeo || size(w, 1) ~= sizeo || size(w, 2) ~= sizeo
            error('Velocity fields do not match the expected sizeo (%d x %d) at time step %d.', sizeo, sizeo, t);
        end

% Compute the velocity at the current position
vectors(:, t) = get_single_velocity(X(:, t), u, v, w, sizeo);

    end
end
