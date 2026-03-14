%% Parameters
dt = 0.05; % Sample time
N = 10; % MPC Horizon
T  = 600; % Simulation steps
nx = 2; % Modified to Polar states [r, phi]
nu = 2; % [v, w]

% Constraints
v_max = 2.0; % macimum linear velocity, m/s
w_max = 1.0; % maximum angular velocity rad/s
a_max = 2.0;     % maximum linear acceleration m/s^2
alpha_max = 4.0; % maximum angular acceleration rad/s^2
dv_max = a_max * dt;
dw_max = alpha_max * dt;
du_max = [dv_max; dw_max]; 
% -----------------------------

% Waypoints extraction from pathIds
num_wp = length(pathIds);
waypoints = zeros(num_wp, 2);
waypoints(1, :) = q_start;

for i = 1:(num_wp - 1)
    curr_node = nodes(pathIds(i)).poly;
    next_node = nodes(pathIds(i+1)).poly;
    intersection_poly = intersect(curr_node, next_node);
    [cx, cy] = centroid(intersection_poly);
    waypoints(i+1, :) = [cx, cy];
end
waypoints(end, :) = q_goal;

% Initialization
x_robot = [q_start(1); q_start(2); 0]; % Robot Cartesian [x; y; theta]
target_idx = 2; 
wp_tol = 0.3;
history_x = zeros(3, T); % Cartesian history for plotting
history_u = zeros(nu, T);
options = optimoptions('quadprog', 'Display', 'off');

%% Simulation Loop
fprintf('Starting MPC Simulation in Polar Coordinates...\n');
for t = 1:T
    % Waypoint Logic
    target = waypoints(target_idx, :)';
    if norm(x_robot(1:2) - target) < wp_tol && target_idx < size(waypoints, 1)
        target_idx = target_idx + 1;
        target = waypoints(target_idx, :)';
    end
    
    % Cartesian to Polar Conversion
    xr = target(1); 
    yr = target(2);
    dx = xr - x_robot(1);
    dy = yr - x_robot(2);
    
    r_state = sqrt(dx^2 + dy^2);
    alpha = atan2(dy, dx);
    phi_state = alpha - x_robot(3);
    phi_state = atan2(sin(phi_state), cos(phi_state)); % Wrap angle to [-pi, pi]
    
    x_polar = [r_state; phi_state];
    
    % Reference velocity for linearization
    if t == 1
        v_ref = 0.1;
    else
        v_ref = max(history_u(1, t-1), 0.1); % Use previous velocity, minimum 0.1
    end
    
    % Jacobian Linearization (Ac and Bc)
    r_safe = max(r_state, 0.05); % Prevent division by zero
    
    Ac = [0,  v_ref * sin(phi_state);
         -v_ref * sin(phi_state) / (r_safe^2),  v_ref * cos(phi_state) / r_safe];
         
    Bc = [-cos(phi_state), 0;
           sin(phi_state) / r_safe, -1];
    
    % Euler Discretization
    A = eye(2) + Ac * dt;
    B = Bc * dt;
    
    % MPC Prediction Matrices
    Phi = zeros(2*N, 2); 
    Gamma = zeros(2*N, 2*N);
    for i = 1:N
        Phi((i-1)*2+1:i*2, :) = A^i;
        for j = 1:i
            Gamma((i-1)*2+1:i*2, (j-1)*2+1:j*2) = (A^(i-j))*B;
        end
    end
    
    % Cost: J = X'*Q*X + U'*R*U
    % Target is to drive r -> 0 and phi -> 0
    Q_b = kron(eye(N), diag([40, 10])); % High penalty on distance, moderate on angle
    R_b = kron(eye(N), diag([0.1, 0.1]));
    
    H = 2 * (Gamma' * Q_b * Gamma + R_b); 
    H = (H + H')/2;
    f = 2 * (x_polar' * Phi' * Q_b * Gamma)';
    
    % Circular Funnel Constraint (r <= r0)
    % Find the current active funnel node
    curr_node_id = pathIds(min(target_idx-1, length(pathIds)));
    
    % Limit max distance based on the circular node's radius to strictly stay inside
    r0 = nodes(curr_node_id).radius * 2.0; % Max possible distance across the funnel
    
    % Build inequality constraints: C_pos * X <= r0
    C_pos = kron(eye(N), [1, 0]); % Extract 'r' states from predictions
    A_funnel = C_pos * Gamma;
    b_funnel = repmat(r0, N, 1) - C_pos * Phi * x_polar;
    
    % MAtrix to get u(k) - u(k-1)
    D = kron(eye(N), eye(2)) - kron(diag(ones(N-1,1), -1), eye(2));

    % Add initial values for u
    if t == 1
        u_prev = [0; 0];
    else
        u_prev = history_u(:, t-1);
    end

    % Cauculate acceleration constraints
    b_upper = repmat(du_max, N, 1);
    b_upper(1:2) = b_upper(1:2) + u_prev;

    b_lower = repmat(du_max, N, 1);
    b_lower(1:2) = b_lower(1:2) - u_prev;

    A_rate = [D; -D];
    b_rate = [b_upper; b_lower];

    % Merge funnel and acceleration constraints
    A_ineq_total = [A_funnel; A_rate];
    b_ineq_total = [b_funnel; b_rate];

    % Solve Quadratic Program
    % lb for velocity is set to 0 to prevent reversing (critical for polar formulation)
    lb = repmat([0; -w_max], N, 1); 
    ub = repmat([v_max;  w_max], N, 1);
    
    [U_opt, ~, exitflag] = quadprog(H, f, A_ineq_total, b_ineq_total, [], [], lb, ub, [], options);
    
    % Apply Control Input
    if exitflag == 1
        u_apply = U_opt(1:2);
    else
        % Proportional fallback in case of optimizer infeasibility near singularity
        u_v = min(v_max, 1.0 * r_state);
        u_w = min(w_max, max(-w_max, 2.0 * phi_state));
        u_apply = [u_v; u_w];
    end
    
    % Update Cartesian Robot State (Kinematics)
    x_robot(1) = x_robot(1) + u_apply(1)*cos(x_robot(3))*dt;
    x_robot(2) = x_robot(2) + u_apply(1)*sin(x_robot(3))*dt;
    x_robot(3) = x_robot(3) + u_apply(2)*dt;
    x_robot(3) = atan2(sin(x_robot(3)), cos(x_robot(3))); % Wrap theta
    
    % Log states and inputs
    history_x(:,t) = x_robot;
    history_u(:,t) = u_apply;
    
    % Early exit condition
    if norm(x_robot(1:2) - waypoints(end,:)') < wp_tol && target_idx == size(waypoints, 1)
        history_x = history_x(:, 1:t);
        history_u = history_u(:, 1:t);
        break;
    end
end

%% Plotting

figure('WindowState','maximized', 'Color','w'); hold on; axis equal;
xlim([W(1) W(2)]); ylim([W(3) W(4)]);
title('Figür 2: Dijkstra - Seçilen Yol ve Kesişim Waypointleri');

for i=1:numel(obs)
    plot(obs{i}, 'FaceColor',[0 0 0], 'FaceAlpha',0.6, 'EdgeColor','none');
end

if ~isempty(pathIds)
    for k = 1:length(pathIds)
        node_idx = pathIds(k);
        plot(nodes(node_idx).poly, 'FaceColor',[1.0 0.85 0.7], 'FaceAlpha',0.40, 'EdgeColor',[1.0 0.5 0.0], 'LineWidth',1.5);
    end

    num_nodes = length(pathIds);
    route_points = zeros(num_nodes + 1, 2); 
    route_points(1, :) = q_start; 

    for k = 1:(num_nodes - 1)
        curr_node = nodes(pathIds(k)).poly;
        next_node = nodes(pathIds(k+1)).poly;

        intersection_poly = intersect(curr_node, next_node);
        [cx, cy] = centroid(intersection_poly);
        route_points(k+1, :) = [cx, cy];
    end

    route_points(end, :) = q_goal; 

    plot(route_points(:,1), route_points(:,2), 'bo', 'MarkerSize',6, 'LineWidth',1.5);
end

plot(q_start(1), q_start(2), 'go', 'MarkerSize',9, 'LineWidth',2);
plot(q_goal(1),  q_goal(2),  'ro', 'MarkerSize',9, 'LineWidth',2);
grid on;

hold on

plot(history_x(1,:), history_x(2,:), 'r-', 'LineWidth', 2);
title('MPC Waypoint Tracking (Fixed)');
legend('Waypoints', 'Robot Path');
grid on; axis equal;
xlabel('X [m]'); ylabel('Y [m]');


t = (0:size(history_u,2)-1) * dt;

figure('Color','w');
subplot(2,1,1);
plot(t, history_u(1,:), 'b', 'LineWidth',1.5);
grid on;
ylabel('u');
title('Control inputs');

subplot(2,1,2);
plot(t, history_u(2,:), 'r', 'LineWidth',1.5);
grid on;
xlabel('Time [s]');
ylabel('w');