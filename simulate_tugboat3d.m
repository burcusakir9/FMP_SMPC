% Simulation script for the 3-DOF Tugboat USV in a straight line
% --- Simulation Parameters ---
T_sim = 100;         % Total simulation time (s)
dt = 0.01;            % Time step (s)
time = 0:dt:T_sim;   % Time vector
% --- Tugboat USV Parameters ---
x_initial = zeros(6, 1); % Initial state: [x y yaw u v r]'
n_desired = [50; 0];    % Desired propeller shaft speeds (rad/s) - adjust for desired speed
% --- Initialize State History ---
x_history = zeros(6, length(time));
x_history(:, 1) = x_initial;
n_history = zeros(2, length(time));
n_history(:, 1) = n_desired;
U_history = zeros(1, length(time));
% --- Simulation Loop ---
x = x_initial;
n = n_desired;
for i = 2:length(time)

    if i < 5000
        n = n_desired;
    else
        n = [0; 50];
    end


    % Calculate the time derivative of the state vector using the 3-DOF model
    xdot = tugboat3d(x, n);
    U = sqrt(x(4)^2 + x(5)^2); % Speed calculation for 3-DOF model
    % Update the state using Euler integration
    x = x + xdot * dt;
    % Store the current state and input
    x_history(:, i) = x;
    n_history(:, i) = n;
    U_history(i) = U;
end
% --- Post-Simulation Analysis and Plotting ---
% Extract states for easier plotting
north = x_history(1, :);
east = x_history(2, :);
yaw = x_history(3, :);
u = x_history(4, :);
v = x_history(5, :);
r = x_history(6, :);
% Plotting the results
figure;
% Position
subplot(2, 2, 1);
plot(east, north);
xlabel('East (m)');
ylabel('North (m)');
title('North-East Position');
grid on;
axis equal;
% Yaw Angle
subplot(2, 2, 2);
plot(time, rad2deg(yaw));
xlabel('Time (s)');
ylabel('Yaw Angle (degrees)');
title('Yaw Angle');
grid on;
% Speed
subplot(2, 2, [3 4]); % Span both bottom subplots
plot(time, x_history(4, :));
xlabel('Time (s)');
ylabel('Speed (m/s)');
title('Speed');
grid on;
sgtitle('3-DOF Tugboat USV Straight Line Simulation');