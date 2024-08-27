function totalEnergy = objective_function(coeffs, t, W, StartLoc, uField, vField, wField, dt)
    try
        % Calculate X
         X = calculateX(t, coeffs, W, StartLoc);
        
        % Check dimensions of X
        disp('Dimensions of X:');
        disp(size(X));
        
        % Get particle locations
        Particles_locations = get_particle_cell_Location(X, 64);
        
        % Check dimensions of Particles_locations
        disp('Dimensions of Particles_locations:');
        disp(size(Particles_locations));
        
        % Get velocity fields at particle locations
        Velocity_fields = get_vectors_at_locations(uField, vField, wField, Particles_locations);
        Velocity_f = Velocity_fields(:,1:end-1);
        
        % Check dimensions of Velocity_fields
        disp('Dimensions of Velocity_fields:');
        disp(size(Velocity_fields));
        disp('Dimensions of Velocity_f:');
        disp(size(Velocity_f));
        
        % Calculate the gap and velocities
        Gap = diff(X, 1, 2);
        Vs = Gap / dt;
        
        % Check dimensions of Gap and Vs
        disp('Dimensions of Gap:');
        disp(size(Gap));
        disp('Dimensions of Vs:');
        disp(size(Vs));
        
        % Calculate the drag forces

        Fdrag = (Velocity_f - Vs) .* vecnorm(Velocity_f - Vs);
        
        % Check dimensions of Fdrag
        disp('Dimensions of Fdrag:');
        disp(size(Fdrag));
        
        % Calculate the total energy
        Totalenergy1 = dot(Fdrag, Gap);
        Totalenergy = sum(Totalenergy1);
        
        totalEnergy = -Totalenergy;
    catch ME
        disp('Error in objective function:');
        disp(ME.message);
        totalEnergy = NaN;
    end
end

