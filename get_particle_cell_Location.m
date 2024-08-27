function cellIndices = get_particle_cell_Location(posArray, numCells)
    % Number of time steps
    numTimeSteps = size(posArray, 2);
    
    % Initialize the array to hold cell indices
    cellIndices = zeros(3, numTimeSteps);
    
    % Iterate through each time step
    for t = 1:numTimeSteps
        % Get the particle position at the current time step
        x = posArray(1, t);
        y = posArray(2, t);
        z = posArray(3, t);
        
        % Convert the continuous position to discrete cell indices
        % Scale the position to the range of the number of cells and add 1
        % because MATLAB indices start at 1
        i = min(max(floor(x * numCells) + 1, 1), numCells);
        j = min(max(floor(y * numCells) + 1, 1), numCells);
        k = min(max(floor(z * numCells) + 1, 1), numCells);
        
        % Store the cell indices
        cellIndices(:, t) = [i; j; k];
    end
end
