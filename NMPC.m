%{
TODO:
- Spiralling around goal point
- Fragile bahıvaior such as meaningless infeasibility in the solution
- Not working with MPC params
%}
%% Parameters

dt = 0.05; % Sample time
N = 5;
T  = 600; % Simnulation steps
nx = 3; % [px, py, theta]
nu = 2; % [v, w]
v_max = 10.0; % m/s
w_max = 1.0; % rad/s

% Waypoints extraction from pathIds
num_wp = length(pathIds);
waypoints = zeros(num_wp, 2);
waypoints(1, :) = q_start;

% Add intersection points to the waypoint list
for i = 1:(num_wp - 1)
    intersection_poly = intersect(nodes(pathIds(i)).poly, nodes(pathIds(i+1)).poly);
    [cx, cy] = centroid(intersection_poly);
    waypoints(i+1, :) = [cx, cy];
end
waypoints(end+1, :) = q_goal; 

% Initialization
x_robot = [q_start(1); q_start(2); 0]; % [x; y; theta]
target_idx = 2; 
wp_tol = 0.5;  % waypoint tolerance to skip to the next
history_x = zeros(nx, T);
history_u = zeros(nu, T);
options = optimoptions('quadprog','Display','off');

%% Simulation Loop

fprintf('Starting NMPC Simulation...\n');

for t = 1:T

     % Waypoint Logic
    target = waypoints(target_idx, :)';
    if norm(x_robot(1:2) - target) < wp_tol && target_idx < size(waypoints, 1)
        target_idx = target_idx + 1;
        target = waypoints(target_idx, :)';
    end

    % Extract Funnel Constraints
    current_node_id = pathIds(min(target_idx-1, length(pathIds)));
    curr_poly = nodes(current_node_id).poly;
    
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
        normal = [edge_vec(2), -edge_vec(1)]; 
        
        if dot(normal, node_centroid - p1) > 0
            normal = -normal; 
        end
        normal = normal / norm(normal);
        
        A_poly(k, :) = normal;
        b_poly(k) = dot(normal, p1);
    end

    % NMPC Setup 
    Q = diag([10, 10, 0.5]); 
    R = diag([0.1, 0.1]);
    
    U0 = zeros(nu * N, 1); 
    
    lb = repmat([-v_max; -w_max], N, 1);
    ub = repmat([ v_max;  w_max], N, 1);
    
    xr = target(1);
    yr = target(2);

    theta_ref = atan2(yr - x_robot(2), xr - x_robot(1));   % fixed reference heading
    
    cost_fun = @(U) nmpc_cost_lq(U, x_robot, xr, yr, theta_ref, Q, R, N, dt);
    nonlcon_fun = @(U) nmpc_constraints(U, x_robot, A_poly, b_poly, N, dt);
    
    opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', 'MaxIterations', 200);
    
    % Solve
    [U_opt, ~, exitflag] = fmincon(cost_fun, U0, [], [], [], [], lb, ub, nonlcon_fun, opts);
    
    if exitflag <= 0
        disp(['NMPC Infeasible at step ' num2str(t) ' - applying brakes.']);
        u_apply = [0;0];
    else
        u_apply = U_opt(1:2);
    end
    
    % Apply to Real Plant 
    x_robot(1) = x_robot(1) + u_apply(1) * cos(x_robot(3)) * dt;
    x_robot(2) = x_robot(2) + u_apply(1) * sin(x_robot(3)) * dt;
    x_robot(3) = x_robot(3) + u_apply(2) * dt;
    x_robot(3) = atan2(sin(x_robot(3)), cos(x_robot(3)));
    
    history_x(:, t) = x_robot;
    history_u(:, t) = u_apply;

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

%% ===================== NMPC FUNCTIONS ===========================
function J = nmpc_cost_lq(U, x0, xr, yr, theta_ref, Q, R, N, dt)
    J = 0;
    x = x0;

    for i = 1:N
        v = U(2*i-1);
        w = U(2*i);

        % Predict nonlinear dynamics
        x(1) = x(1) + v * cos(x(3)) * dt;
        x(2) = x(2) + v * sin(x(3)) * dt;
        x(3) = x(3) + w * dt;
        x(3) = atan2(sin(x(3)), cos(x(3)));

        % Target-relative error state (same idea as your linear MPC)
        e1 = x(1) - xr;
        e2 = x(2) - yr;
        e3 = x(3) - theta_ref;
        e3 = atan2(sin(e3), cos(e3));     % wrap angle error

        e = [e1; e2; e3];

        % Quadratic stage cost
        J = J + e' * Q * e + [v; w]' * R * [v; w];
    end
end

function [c, ceq] = nmpc_constraints(U, x0, A_poly, b_poly, N, dt)
    ceq = []; 
    c = zeros(size(A_poly,1) * N, 1);
    x = x0;
    
    for i = 1:N
        v = U(2*i-1);
        w = U(2*i);
        
        x(1) = x(1) + v * cos(x(3)) * dt;
        x(2) = x(2) + v * sin(x(3)) * dt;
        x(3) = x(3) + w * dt;
        x(3) = atan2(sin(x(3)), cos(x(3))); 
        
        idx_start = (i-1)*size(A_poly,1) + 1;
        idx_end = i*size(A_poly,1);
        
        c(idx_start:idx_end) = A_poly * x(1:2) - b_poly;
    end
end