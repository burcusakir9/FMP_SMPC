%% Parameters

dt_sim = 0.01; % Sample time of simulation
dt_mpc = 0.05; % Sample time of MPC
N = 10; % MPC Horizon
sim_time = 50; % seconds
T  = sim_time / dt_sim; % Simulation steps
mpc_interval = round(dt_mpc / dt_sim);

nx = 2; % Modified to Polar states [r, phi]
nu = 2; % [v, w]

% Constraints
v_max = 2.0; % macimum linear velocity, m/s
w_max = 1.0; % maximum angular velocity rad/s
a_max = 2.0;     % maximum linear acceleration m/s^2
alpha_max = 4.0; % maximum angular acceleration rad/s^2
dv_max = a_max * dt_mpc;
dw_max = alpha_max * dt_mpc;
du_max = [dv_max; dw_max]; 

wp_tol = 0.3;
goal_point_tol =  0.01;

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
x_robot = [q_start(1); q_start(2); 0]; % Robot cartesian [x; y; theta]
target_idx = 2; 
history_x = zeros(3, T); % Cartesian history for plotting
history_u = zeros(nu, T);

options = optimoptions('quadprog', 'Display', 'off');

%% Simulation Loop

fprintf('Starting MPC Simulation in Polar Coordinates...\n');

for t = 1:T
    if mod(t-1, mpc_interval) == 0

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
        A = eye(2) + Ac * dt_mpc;
        B = Bc * dt_mpc;
        
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
        Q_b = kron(eye(N), diag([40, 10])); % High penalty on distance, moderate on angle
        R_b = kron(eye(N), diag([0.1, 0.1]));
        
        H = 2 * (Gamma' * Q_b * Gamma + R_b); 
        H = (H + H')/2;
        f = 2 * (x_polar' * Phi' * Q_b * Gamma)';
        
        % Active nodes
        active_nodes = nodes(pathIds);
        current_poly_idx = find_node_for_point([xr, yr], active_nodes); 
        curr_node = active_nodes(current_poly_idx);
        
        % target to center distance
        d_target_to_center = norm([xr, yr] - curr_node.c);
        
        r0 = max(curr_node.radius - d_target_to_center, 0.05); % Avoid division by zero
    
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
    end
    
    % Update Cartesian Robot State (Kinematics)
    x_robot(1) = x_robot(1) + u_apply(1)*cos(x_robot(3))*dt_sim;
    x_robot(2) = x_robot(2) + u_apply(1)*sin(x_robot(3))*dt_sim;
    x_robot(3) = x_robot(3) + u_apply(2)*dt_sim;
    x_robot(3) = atan2(sin(x_robot(3)), cos(x_robot(3))); % Wrap theta
    
    % Log states and inputs
    history_x(:,t) = x_robot;
    history_u(:,t) = u_apply;
    
    % Early exit condition
    if norm(x_robot(1:2) - q_goal') < goal_point_tol
        history_x = history_x(:, 1:t);
        history_u = history_u(:, 1:t);
        break;
    end
end
%% Helper functions
function idx = find_node_for_point(pt, nodes)
    min_dist = inf;
    closest_idx = 1;
    
    for i = 1:numel(nodes)
        % Check if the point is strictly inside the polygon
        if isinterior(nodes(i).poly, pt(1), pt(2))
            idx = i;
            return;
        end
        
        % Distance fallback calculation in case of spline overshoot
        dist = norm(pt - nodes(i).c);
        if dist < min_dist
            min_dist = dist;
            closest_idx = i;
        end
    end
    
    % Return the closest node if the point is outside all polygons
    idx = closest_idx;
    
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
title('MPC Waypoint Tracking');
legend('Waypoints', 'Robot Path');
grid on; axis equal;
xlabel('X [m]'); ylabel('Y [m]');


t = (0:size(history_u,2)-1) * dt_sim;

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