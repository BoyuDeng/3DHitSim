function [uField, vField, wField] = ChangeU(u, v, w, fac)



% Preallocate cells to store the multiplied data
    uField = cell(size(u));
    vField= cell(size(v));
    wField = cell(size(w));


for t = 1:length(u)
    uMat = cell2mat(u{t}) .* fac;
    vMat = cell2mat(v{t}) .* fac;
    wMat = cell2mat(w{t}) .* fac;

    if isnumeric(uMat) && isnumeric(vMat) && isnumeric(wMat)
        uField{t} = mat2cell(uMat, size(uMat,1), size(uMat,2), size(uMat,3));
        vField{t} = mat2cell(vMat, size(vMat,1), size(vMat,2), size(vMat,3));
        wField{t} = mat2cell(wMat, size(wMat,1), size(wMat,2), size(wMat,3));
    else
        error('Each cell must contain numeric data.');
    end
end
