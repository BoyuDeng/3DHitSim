for i = 1:length(G)
    X = X14unlim(t, resultsU8_8{i}.coeffs, resultsU8_8{i}.W);
    Zdiff(i) = X(3, end)-X(3, 1);
end