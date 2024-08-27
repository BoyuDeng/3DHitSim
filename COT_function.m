function COT = COT_function(coeffs, t, W, StartLoc, uField, vField, wField, dt, p,U)
    try


        tau_p = 1.22;
        G = 9.8*tau_p/U;
        St = 0;
        Du = calculateDu(t,coeffs);
        % Calculate X
         X = calculateX(t, coeffs, W, StartLoc);
        
        % Check dimensions of X
        % disp('Dimensions of X:');
        % disp(size(X));
        
        % Get particle locations
        Particles_locations = get_particle_cell_Location(X, 64);
        
        % Check dimensions of Particles_locations
        % disp('Dimensions of Particles_locations:');
        % disp(size(Particles_locations));
        % 
        % Get velocity fields at particle locations
        Velocity_fields = get_vectors_at_locations(uField, vField, wField, Particles_locations);
        
        % Check dimensions of Velocity_fields
        % disp('Dimensions of Velocity_fields:');
        % disp(size(Velocity_fields));
        % 
        % Calculate the gap and velocities
        Vs = calculateV(t, coeffs, W);
        
        % Check dimensions of Gap and Vs
        % disp(size(Vs));
        % 
        % Calculate the drag forces

        Fdrag =(1/G)* (St*Du-(Velocity_fields - Vs) .* vecnorm(Velocity_fields - Vs).^(p-1))+1;
        
        % Check dimensions of Fdrag
        % disp('Dimensions of Fdrag:');
        % disp(size(Fdrag));
        
        % Calculate the total energy
        
        COT = (G/(W*t(end)))*sum((vecnorm(Fdrag).^2).^(3/4))*dt;
    catch ME
        disp('Error in objective function:');
        disp(ME.message);
        COT = NaN;
    end
end
