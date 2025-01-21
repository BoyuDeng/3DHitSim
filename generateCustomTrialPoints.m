% Custom trial point generator function
function points = generateCustomTrialPoints(numPoints, W)
    % Generate trial points with only the first four coefficients varying
    points = zeros(numPoints, 27);  % 26 coefficients + W
    points(:, 1:4) = -1 + 2 * rand(numPoints, 4);  % First 4 coefficients in [-1, 1]
    points(:, 27) = W;  % Fix W to its initial value
end