function [uField, vField, wField] = CelltoMat(uField, vField, wField)    

% Validate input arguments
    if ~iscell(uField) || ~iscell(vField) || ~iscell(wField)
        error('Input arguments must be cell arrays.');
    end

    % Get the dimensions of the cell arrays
    numCells = numel(uField);
    
    % Assume each cell contains matrices of the same size
    field = uField{1};
    field = cell2mat(field);
    [rows, cols, pa] = size(field);
    
    % Preallocate matrices to store the data
    uFieldMat = zeros(rows, cols, pa, numCells);
    vFieldMat = zeros(rows, cols, pa, numCells);
    wFieldMat = zeros(rows, cols, pa, numCells);

    % Convert cell arrays to matrices using a loop
    for t = 1:numCells
        u = cell2mat(uField{t});
        v = cell2mat(vField{t});
        w = cell2mat(wField{t});

        if isnumeric(u) && isnumeric(v) && isnumeric(w)
            uFieldMat(:,:,:,t) = u;
            vFieldMat(:,:,:,t) = v;
            wFieldMat(:,:,:,t) = w;
        else
            error('Each cell must contain numeric data.');
        end
    end

    % Calculate the mean velocities across all dimensions


    % Calculate the velocity fluctuations
    uField = uFieldMat;
    vField = vFieldMat;
    wField = wFieldMat;