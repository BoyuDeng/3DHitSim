function X = QFode45(uField, vField, wField, t, dt, sizeo, z, x)
    X0 = [0.5+x, 0, 0.5+z];  % Set initial position
    dt = 2*dt;
    t = 0:dt:t(end);
    % Call ode45
    n = length(t);
    X = zeros(length(X0), n);
    X(:,1) = X0;
  
    for i = 1:n-1
        k1 = get_time_vel(uField, vField, wField, sizeo,i, X(:,i));
        k2 = get_time_vel(uField, vField, wField, sizeo,(i + 1), X(:,i) + dt*k1/2);
        k3 = get_time_vel(uField, vField, wField, sizeo,(i + 1), X(:,i) + dt*k2/2);
        k4 = get_time_vel(uField, vField, wField, sizeo, (i + 2), X(:,i) + dt*k3);
        X(:,i+1) = X(:,i) + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
    end


    % Plot pathline
figure;
plot3(X(1,:)/0.17, X(2,:)/0.17, X(3,:)/0.17, 'b', 'LineWidth', 1.5); % Blue line, no markers
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Particle Pathline');
grid on;
axis equal;  % Keep aspect ratio consistent
view(3);  % Ensure a good 3D perspective