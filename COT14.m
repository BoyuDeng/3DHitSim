function [COT ,Ft, G, alitotal, Dterm] = COT14(coeffs, t, W, uField, vField, wField, dt, p,U, B, G)
    try


       
        %St = tau_p*U/0.2;
        St = 0;
        Du = 0;   %Du14(t,coeffs);
        % Calculate X
        X = X14unlim(t, coeffs, W);
        
        
 
        Velocity_fields = get_vel(uField,vField,wField,64,X);
        

        Vs = V14unlim(t, coeffs, W);
        

        Ft =(1/G)* (St*Du-(Velocity_fields./U - Vs./U) .* vecnorm(Velocity_fields./U - Vs./U).^(p-1));
        ali = (Velocity_fields./U - Vs./U);
        Ft(3,:) = Ft(3,:) + B;
        alitotal = vecnorm(ali);
        Dterm = (Velocity_fields./U - Vs./U);
        
        COT = ((G*U)/(X(2,end)))*sum((vecnorm(Ft).^2).^(3/4))*(dt);


if p == 1
    COT = abs(COT / 1.611);
elseif p == 2
    Eqf = (G / sqrt(G * sqrt(1/2))) * (1 + sqrt(1/2))^(3/4);
    COT = abs(COT / Eqf);
end

    catch ME
        disp('Error in objective function:');
        disp(ME.message);
        COT = NaN;
    end
end
