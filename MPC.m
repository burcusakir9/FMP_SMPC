%{
TODO:
- Spiralling around goal point
- Chattering in the input
%}

%% Parameters

dt = 0.05; % Sample time
N = 10; % MPC Horizon
T  = 600; % Simnulation steps
nx = 3; % [px, py, theta]
nu = 2; % [v, w]
v_max = 2.0; % m/s
w_max = 1.0; % rad/s

% Waypoints extraction from pathIds
num_wp = length(pathIds);
waypoints = zeros(num_wp, 2);
waypoints(1, :) = q_start;

for i = 1:(num_wp - 1)
    intersection_poly = intersect(nodes(pathIds(i)).poly, nodes(pathIds(i+1)).poly);
    [cx, cy] = centroid(intersection_poly);
    waypoints(i+1, :) = [cx, cy];
end
waypoints(end, :) = q_goal;


% Initialization
x_robot = [q_start(1); q_start(2); 0]; % [x; y; theta]
target_idx = 2; 
wp_tol = 0.3;
history_x = zeros(nx, T);
history_u = zeros(nu, T);
options = optimoptions('quadprog', 'Display', 'off');


%% Simulation Loop

fprintf('Starting MPC Simulation...\n');

for t = 1:T
    % Waypoint Logic
    target = waypoints(target_idx, :)';
    if norm(x_robot(1:2) - target) < wp_tol && target_idx < size(waypoints, 1)
        target_idx = target_idx + 1;
        target = waypoints(target_idx, :)';
    end
    
    % Linearization
    yaw = x_robot(3);
    if t == 1
        v_ref = 0.1;
    else
        v_ref = max(abs(history_u(1, t-1)), 0.1);
    end
    
    A = [1, 0, -v_ref*sin(yaw)*dt; 
         0, 1,  v_ref*cos(yaw)*dt; 
         0, 0,  1];
    B = [cos(yaw)*dt, 0; 
         sin(yaw)*dt, 0; 
         0,           dt];
    
    % Shift to target-relative frame
    xr = target(1); yr = target(2);
    x_bar = [x_robot(1)-xr; x_robot(2)-yr; 0];
    target_theta = atan2(yr - x_robot(2), xr - x_robot(1));
    d_th = target_theta - x_robot(3);
    x_bar(3) = -atan2(sin(d_th), cos(d_th));
    
    % MPC Prediction Matrices
    Phi = zeros(3*N, 3); 
    Gamma = zeros(3*N, 2*N);
    for i = 1:N
        Phi((i-1)*3+1:i*3, :) = A^i;
        for j = 1:i
            Gamma((i-1)*3+1:i*3, (j-1)*2+1:j*2) = (A^(i-j))*B;
        end
    end
    
    % Cost: J = X'*Q*X + U'*R*U
    Q_b = kron(eye(N), diag([40, 40, 1])); 
    R_b = kron(eye(N), diag([0.1, 0.1]));
    H = 2 * (Gamma' * Q_b * Gamma + R_b); H = (H + H')/2;
    f = 2 * (x_bar' * Phi' * Q_b * Gamma)';
    
    % Funnel Constraints (Ax <= b)
    curr_node_id = pathIds(min(target_idx-1, length(pathIds)));
    [vx, vy] = boundary(nodes(curr_node_id).poly);
    valid = ~isnan(vx); vx=vx(valid); vy=vy(valid);
    if ~isempty(vx) && vx(1)==vx(end), vx(end)=[]; vy(end)=[]; end
    
    num_e = length(vx); A_p = zeros(num_e, 2); b_p = zeros(num_e, 1);
    center = nodes(curr_node_id).c;
    for k = 1:num_e
        p1 = [vx(k), vy(k)]; p2 = [vx(mod(k,num_e)+1), vy(mod(k,num_e)+1)];
        nv = [p2(2)-p1(2), -(p2(1)-p1(1))];
        if dot(nv, center - p1) < 0, nv = -nv; end
        nv = nv/norm(nv);
        A_p(k,:) = -nv; b_p(k) = dot(-nv, p1);
    end
    
    C_pos = kron(eye(N), [1 0 0; 0 1 0]);
    A_ineq = kron(eye(N), A_p) * C_pos * Gamma;
    b_ineq = repmat(b_p, N, 1) - kron(eye(N), A_p) * (C_pos*Phi*x_bar + repmat([xr;yr], N, 1));
    
    % Solve
    lb = repmat([-v_max; -w_max], N, 1); 
    ub = repmat([ v_max;  w_max], N, 1);
    [U_opt, ~, exitflag] = quadprog(H, f, A_ineq, b_ineq, [], [], lb, ub, [], options);
    
    % Apply
    if exitflag == 1, u_apply = U_opt(1:2); else, u_apply = [0.1; 0]; end
    
    x_robot(1) = x_robot(1) + u_apply(1)*cos(x_robot(3))*dt;
    x_robot(2) = x_robot(2) + u_apply(1)*sin(x_robot(3))*dt;
    x_robot(3) = x_robot(3) + u_apply(2)*dt;
    x_robot(3) = atan2(sin(x_robot(3)), cos(x_robot(3)));
    
    history_x(:,t) = x_robot;
    history_u(:,t) = u_apply;
    
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