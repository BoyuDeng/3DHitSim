function [COT,Fdrag, G, ali] = COT14(coeffs, t, W, uField, vField, wField, dt, p,U, B)
    try


        tau_p = 1;
        G = 9.8*tau_p/U;
        %St = tau_p*U/0.2;
        St = 0;
        Du = Du14(t,coeffs);
        % Calculate X
        X = X14(t, coeffs, W);
        
        % Check dimensions of X
        % disp('Dimensions of X:');
        % disp(size(X));
        
        % Get particle locations
        %Particles_locations = get_particle_cell_Location(X, 64);
        
        % Check dimensions of Particles_locations
        % disp('Dimensions of Particles_locations:');
        % disp(size(Particles_locations));
        % 
        % Get velocity fields at particle locations
        %Velocity_fields = get_vectors_at_locations(uField, vField, wField, Particles_locations);
        Velocity_fields = get_vel(uField,vField,wField,64,X);
        
        % Check dimensions of Velocity_fields
        % disp('Dimensions of Velocity_fields:');
        % disp(size(Velocity_fields));
        % 
        % Calculate the gap and velocities
        Vs = V14(t, coeffs, W);
        
        % Check dimensions of Gap and Vs
        % disp(size(Vs));
        % 
        % Calculate the drag forces

        Fdrag =(1/G)* (St*Du-(Velocity_fields./U - Vs./U) .* vecnorm(Velocity_fields./U - Vs./U).^(p-1));
        ali = abs(Velocity_fields./U - Vs./U);
        Fdrag(3,:) = Fdrag(3,:) - B;
        
        COT = (G/(W*t(end)))*sum((vecnorm(Fdrag).^2).^(3/4))*(dt/t(end));

        
    catch ME
        disp('Error in objective function:');
        disp(ME.message);
        COT = NaN;
    end
end
