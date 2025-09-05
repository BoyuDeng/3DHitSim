%% PARAMETERS -----------------------------------------------------------
L          = 0.20;                       % normalisation length
trajIdx    = [1 2 3];                 % <-- choose which “X” data sets you want
numTraj    = numel(trajIdx);             %   to overlay (3 in this case)
colors     = lines(numTraj);             % distinct colours for each

%% PRE-ALLOCATION --------------------------------------------------------
num_time_steps = length(t);
X   = zeros(3,num_time_steps);           % workspace for X14unlim output
U_  = zeros(1,num_time_steps);
V_  = zeros(1,num_time_steps);
W_  = zeros(1,num_time_steps);

%% FIGURE ---------------------------------------------------------------
figure
hold on; grid on; axis equal
title('Trajectories G=0.2')
xlabel('X/L'); ylabel('Y/L'); zlabel('Z/L')

% LOOP OVER TRAJECTORIES ----------------------------------------------
for kk = 1:numTraj
     dataIdx          = trajIdx(kk);   
    % Define the file path dynamically
    
    
    % ----- 1.  build the trajectory ------------------------------------
                            % pick the data set
    coeffStruct      = resultsallU8_2{dataIdx}{23}.coeffs;
    wStruct          = resultsallU8_2{dataIdx}{23}.W;
    X(:,:)           = X14unlim(t, coeffStruct, wStruct);      % [3 × Nt] array
    
    X_k = X(1,:)./L;                                           % normalised coords
    Y_k = X(2,:)./L;
    Z_k = X(3,:)./L;
    
    plot3(X_k,Y_k,Z_k,'LineWidth',2,'Color',colors(kk,:), ...
          'DisplayName',['Trajectory ' num2str(kk)]);
    
    % ----- 2.  velocity field arrows -----------------------------------
    V          = get_vel(uField,vField,wField,64,X);           % [3 × Nt]
    V_scale    = 1;                                            % arrow scaling
    U_(:)      =  V(1,:).*V_scale./L;
    V_(:)      =  V(2,:).*V_scale./L;
    W_(:)      =  V(3,:).*V_scale./L;
    
    step = 10;                                                 % draw every 10th
    %quiver3(X_k(1:step:end),Y_k(1:step:end),Z_k(1:step:end), ...
            %
end

% AXES & LEGEND --------------------------------------------------------
xlim([0 1/L]); ylim([0 1/L]); zlim([0 1/L]);
view(3)
legend('show','Location','best')
hold off
