parfor i = 1:N
    [results2Z4_5{i}.coeffs, results2Z4_5{i}.W, results2Z4_5{i}.energy, fvals2Z4_5{i}] = ...
        optimization27(t, U, uField, vField, wField, dt, p, U, Forcing, 1000,5, G(i),0,20,1,1,0);
end 

parfor i = 1:N
    [results2Z4_10{i}.coeffs, results2Z4_10{i}.W, results2Z4_10{i}.energy, fvals2Z4_10{i}] = ...
        optimization27(t, results2Z4_5{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000,10, G(i),0,20,1,1,results2Z4_5{i}.coeffs);
end

parfor i = 1:N
    [results2Z4_15{i}.coeffs, results2Z4_15{i}.W, results2Z4_15{i}.energy, fvals2Z4_15{i}] = ...
        optimization27(t, results2Z4_10{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000,15, G(i),0, 20,1,1,results2Z4_10{i}.coeffs);
end

parfor i = 1:N
    [results2Z4_20{i}.coeffs, results2Z4_20{i}.W, results2Z4_20{i}.energy, fvals2Z4_20{i}] = ...
        optimization27(t, results2Z4_15{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000,20, G(i),0, 20,1,1,results2Z4_15{i}.coeffs);
end

parfor i = 1:N
    [results2Z4_25{i}.coeffs, results2Z4_25{i}.W, results2Z4_25{i}.energy, fvals2Z4_25{i}] = ...
        optimization27(t, results2Z4_20{i}.W, uField, vField, wField, dt, p, U, Forcing, 1000,25, G(i),0, 20,1,1,results2Z4_20{i}.coeffs);
end

%%
for i = 1:14
coeffs = results2Z4_10{i}.coeffs;
W = results2Z4_10{i}.W;


        X = X14unlim(t, coeffs, W);

        Velocity_fields = get_vel(uField,vField,wField,64,X);
        
        Vs = V14unlim(t, coeffs, W);
        

end

VV = abs(Velocity_fields);
mean(VV(2,:))