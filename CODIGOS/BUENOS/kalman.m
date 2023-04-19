% RADAR TRACKING. KALMAN FILTERING

% System parameters
T = 0.1;
Q = [T^5/20 T^4/8 T^3/6 ; T^4/8 T^3/3 T^2/2 ; T^3/6 T^2/2 T];
F = [1 T T^2/2 ; 0 1 T ; 0 0 1];
H = [1 0 0 ; 0 1 0];
R = eye(2)*0.1;

steps = 100;
z = zeros(2,steps);

% Intialization
x = [z(1,1) ; z(2,1) ; 0];

% TODO Determinar valor sigmas
P = eye(3)*1;

stored_measurements = zeros(3,step);
measurement_region = zeros(1,step);

for k = 1:steps
    % Prediction
    x = F*x;
    P = F*P*F' + Q;
    v = z(:,k) - H*x;
    
    % Update
    S = H*P*H' + R;
    gate = v'*(1/S)*v;
    K = (P*H')/S;
    x = x + K*v;
    P = P - K*H*P;
    
    % Store results
    stored_measurements(:,k) = x;
    measurement_region(k) = gate;
end




