%{
TODO:
- 
%}
%% Parameters

dt = 0.05; % sample time

N = 5;

T  = 300; % Simnulation steps
nx = 3; % [px, py, theta]
nu = 2; % [v, w]

% Input limits
v_max = 10.0; % m/s
w_max = 5.0; % rad/s

% Waypoints extraction from pathIds
num_wp = length(pathIds);
waypoints = zeros(num_wp, 2);

% Add start point to the waypoint list
waypoints(1, :) = q_start;

% Add intersection points to the waypoint list
for i = 1:(num_wp - 1)
    curr_node = nodes(pathIds(i)).poly;
    next_node = nodes(pathIds(i+1)).poly;
    
    % Hedef = İki düğümün kesişiminin merkezi (Her iki poligonun da içindedir!)
    intersection_poly = intersect(curr_node, next_node);
    [cx, cy] = centroid(intersection_poly);
    waypoints(i+1, :) = [cx, cy];
end

% Add end point to the waypoint list
waypoints(end+1, :) = q_goal; 

wp_tol = 0.5;  % waypoint tolerance to skip to the next
current_wp_idx = 2; 
target = waypoints(current_wp_idx, :)'; 
x_true = [waypoints(1,1); waypoints(1,2); 0]; % Robot starts at the start point

% Storage
history_x = zeros(nx, T);
history_u = zeros(nu, T);

options = optimoptions('quadprog','Display','off');

%% Simulation Loop
for t = 1:T
    % Waypoint Logic
    dist_to_target = norm(x_true(1:2) - target(1:2));

    is_last_waypoint = (current_wp_idx == size(waypoints, 1));
    
    if dist_to_target < wp_tol
        if is_last_waypoint
            disp(['Final destination reached. Station Keeping... (Step: ' num2str(t) ')']);
            
            % Stop
            if dist_to_target < 0.1 && norm(history_u(:, t-1)) < 0.05
                 disp('Robot has fully stopped. Simulation End.');
                 break;
            end
        else
            % If not last waypoint, iterate
            disp(['Waypoint ' num2str(current_wp_idx) ' reached!']);
            current_wp_idx = current_wp_idx + 1;
            target = waypoints(current_wp_idx, :)';
        end
    end

    % Linearization
    th = x_true(3);
    v_ref = 1.0; 
    
    A = [1 0 -v_ref*sin(th)*dt;
         0 1  v_ref*cos(th)*dt;
         0 0  1];
    B = [cos(th)*dt  0;
         sin(th)*dt  0;
         0           dt];

   % Local Coordinate Shift
    xr = target(1);
    yr = target(2);
    
    dx = xr - x_true(1);
    dy = yr - x_true(2);
    raw_target_theta = atan2(dy, dx);
    
    % Local state
    x_bar = zeros(3, 1);
    x_bar(1) = x_true(1) - xr;
    x_bar(2) = x_true(2) - yr;
    
    if dist_to_target < 0.5
        x_bar(3) = 0; 
    else
        delta_theta = raw_target_theta - x_true(3);
        x_bar(3) = -atan2(sin(delta_theta), cos(delta_theta)); 
    end
    
    % MPC Matrices
    Q = diag([10, 10, 0.5]); 
    R = diag([0.1, 0.1]);
    
    Phi = zeros(nx*N, nx);
    Gamma = zeros(nx*N, nu*N);
    temp_A = A;
    for i = 1:N
        Phi((i-1)*nx+1:i*nx, :) = temp_A;
        for j = 1:i
             if j == i, mult = B; else, mult = (A^(i-j))*B; end
             Gamma((i-1)*nx+1:i*nx, (j-1)*nu+1:j*nu) = mult;
        end
        temp_A = temp_A * A;
    end
    
    Q_bar = kron(eye(N), Q);
    R_bar = kron(eye(N), R);
    
    % Origin in local reference frame
    Ref = zeros(nx*N, 1);
    
    H_qp = Gamma' * Q_bar * Gamma + R_bar;
    f_qp = (x_bar' * Phi' * Q_bar * Gamma - Ref' * Q_bar * Gamma)';
    
    % Funnel Constraints
    current_node_id = pathIds(current_wp_idx - 1);
    curr_poly = nodes(current_node_id).poly;
    
    % Poligon limits [A_poly * x <= b_poly] 
    [vx, vy] = boundary(curr_poly);
    if isnan(vx(end)), vx(end)=[]; vy(end)=[]; end
    if vx(1)==vx(end) && vy(1)==vy(end), vx(end)=[]; vy(end)=[]; end
    
    num_edges = length(vx);
    A_poly = zeros(num_edges, 2);
    b_poly = zeros(num_edges, 1);
    node_centroid = nodes(current_node_id).c;
    
    for k = 1:num_edges
        p1 = [vx(k), vy(k)];
        p2 = [vx(mod(k, num_edges)+1), vy(mod(k, num_edges)+1)];
        edge_vec = p2 - p1;
        normal = [edge_vec(2), -edge_vec(1)]; % Normal vektör
        
        if dot(normal, node_centroid - p1) > 0
            normal = -normal; 
        end
        normal = normal / norm(normal);
        
        A_poly(k, :) = normal;
        b_poly(k) = dot(normal, p1);
    end
    
    % Local frame shift: A_poly * (x_bar + x_r) <= b_poly => A_poly * x_bar <= b_poly - A_poly * x_r
    A_rep = kron(eye(N), A_poly);
    b_rep = repmat(b_poly - A_poly * [xr; yr], N, 1);
    
    C_pos = kron(eye(N), [1 0 0; 0 1 0]); 
    
    % Inequality constraints: A_ineq * U <= b_ineq
    A_ineq = A_rep * C_pos * Gamma;
    b_ineq = b_rep - A_rep * C_pos * Phi * x_bar;
    
    % Input constraints
    lb = repmat([-v_max; -w_max], N, 1);
    ub = repmat([ v_max;  w_max], N, 1);
    
    % Solve
    [U_opt, ~, exitflag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, [], [], lb, ub, [], options);
    
    if exitflag ~= 1
        disp(['QP Infeasible at step ' num2str(t) ' - applying brakes.']);
        u_apply = [0;0];
    else
        u_apply = U_opt(1:2);
    end
    
    % Apply to Real Plant
    x_true(1) = x_true(1) + u_apply(1) * cos(x_true(3)) * dt;
    x_true(2) = x_true(2) + u_apply(1) * sin(x_true(3)) * dt;
    x_true(3) = x_true(3) + u_apply(2) * dt;
    x_true(3) = atan2(sin(x_true(3)), cos(x_true(3)));
    
    history_x(:, t) = x_true;
    history_u(:, t) = u_apply;

end

% Truncate
if t < T
    final_step = t - 1;
else
    final_step = T;
end

history_x = history_x(:, 1:final_step);
history_u = history_u(:, 1:final_step);

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
% plot(waypoints(:,1), waypoints(:,2), 'k--s', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
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