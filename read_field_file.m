function [dims, valnames, values, varnames, variables] = read_field_file(filename)
    % Open the file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file %s', filename);
    end

    % Read dimensions
    nx = fread(fid, 1, 'int32');
    ny = fread(fid, 1, 'int32');
    nz = fread(fid, 1, 'int32');
    nval = fread(fid, 1, 'int32');
    nvar = fread(fid, 1, 'int32');
    
    % Store dimensions in a struct
    dims = struct('nx', nx, 'ny', ny, 'nz', nz, 'nval', nval, 'nvar', nvar);

    % Display dimensions
    fprintf('Dimensions: nx = %d, ny = %d, nz = %d\n', nx, ny, nz);
    fprintf('Number of values: %d\n', nval);
    fprintf('Number of variables: %d\n', nvar);

    % Read value names
    valnames = strings(nval, 1);
    for i = 1:nval
        valnames(i) = strtrim(fread(fid, 8, '*char')');
    end
    fprintf('Value names: %s\n', strjoin(valnames, ', '));

    % Read values
    values = fread(fid, nval, 'double');
    disp('Values:');
    disp(values');

    % Read variable names
    varnames = strings(nvar, 1);
    for i = 1:nvar
        varnames(i) = strtrim(fread(fid, 8, '*char')');
    end
    fprintf('Variable names: %s\n', strjoin(varnames, ', '));

    % Read variables
    variables = cell(nvar, 1);
    expected_size = nx * ny * nz;
    for i = 1:nvar
        data = fread(fid, expected_size, 'double');
        if length(data) ~= expected_size
            warning('Variable %d expected %d elements, but got %d elements.', i, expected_size, length(data));
            data = [data; zeros(expected_size - length(data), 1)]; % Pad with zeros if needed
        end
        variables{i} = reshape(data, [nx, ny, nz]);
    end

    % Close the file
    fclose(fid);
end
