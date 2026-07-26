%% Parameters

dt_sim = 0.01;          % Simulation sample time
dt_mpc = 0.05;          % MPC sample time
N = 10;                 % MPC horizon
sim_time = 200;         % seconds
T = round(sim_time / dt_sim);
mpc_interval = round(dt_mpc / dt_sim);

nx = 6;                 % [x, y, yaw, u, v, r]
nu = 2;                 % [n1, n2]

% Propeller limits
n_max = 100;
n_min = -100;

% Propeller slew-rate limits (per MPC step)
dn_max = 100;
du_max = [dn_max; dn_max];

% Waypoint / stopping tolerances
wp_tol = 5.0;
goal_point_tol = 2.0;

% MPC weights
Q = diag([10, 10, 200, 5, 5, 50]);
R = diag([0.5, 0.5]);

%% Waypoints extraction from pathIds
num_wp = length(pathIds);
waypoints = zeros(num_wp + 1, 2);

waypoints(1, :) = q_start;

for i = 1:(num_wp - 1)
    curr_node = nodes(pathIds(i)).poly;
    next_node = nodes(pathIds(i+1)).poly;
    intersection_poly = intersect(curr_node, next_node);

    if isempty(intersection_poly) || intersection_poly.NumRegions == 0 || area(intersection_poly) <= 0
        cp = 0.5 * (nodes(pathIds(i)).c + nodes(pathIds(i+1)).c);
    else
        [cx, cy] = centroid(intersection_poly);
        cp = [cx, cy];
    end

    waypoints(i+1, :) = cp;
end

waypoints(end, :) = q_goal;

%% Initialization
x_robot = [q_start(1); q_start(2); 0; 0; 0; 0];

if size(waypoints,1) >= 2
    target_idx = 2;     % start tracking first real waypoint after q_start
else
    target_idx = 1;
end

history_x = zeros(nx, T);
history_u = zeros(nu, T);

u_prev = [0; 0];
u_apply = [0; 0];

qp_options = optimoptions('quadprog', 'Display', 'off');
nlp_options = optimoptions('fmincon', ...
    'Display', 'off', ...
    'Algorithm', 'sqp', ...
    'MaxIterations', 100, ...
    'StepTolerance', 1e-8, ...
    'ConstraintTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6);

fprintf('Starting tugboat pure linear deviation-form MPC waypoint tracking with direct circular funnel constraints.\n');

%% Simulation Loop
for t = 1:T

    if mod(t-1, mpc_interval) == 0

        %% Waypoint update
        target = waypoints(target_idx, :)';

        if norm(x_robot(1:2) - target) < wp_tol && target_idx < size(waypoints, 1)
            target_idx = target_idx + 1;
            target = waypoints(target_idx, :)';
        end

        %% Reference state
        dx = target(1) - x_robot(1);
        dy = target(2) - x_robot(2);

        yaw_ref = atan2(dy, dx);
        yaw_err = atan2(sin(yaw_ref - x_robot(3)), cos(yaw_ref - x_robot(3)));
        yaw_ref = x_robot(3) + yaw_err;

        x_ref = [target(1);
                 target(2);
                 yaw_ref;
                 0;
                 0;
                 0];

        %% Moving operating point for linearization
        x_lin = x_robot;
        u_lin = u_prev;

        if norm(u_lin) < 1e-6
            u_lin = [20; 20];
        end

        %% Continuous-time local linearization
        [A_c, B_c] = tugboat3d_linearized(x_lin, u_lin);

        %% Discretization (forward Euler, pure linear deviation model)
        % delta_x(k+1) = A*delta_x(k) + B*delta_u(k)
        A = eye(nx) + A_c * dt_mpc;
        B = B_c * dt_mpc;

        %% Deviation initial condition and reference
        delta_x0    = x_robot - x_lin;   % zero by construction
        delta_x_ref = x_ref   - x_lin;

        %% Prediction matrices for linear model
        Phi   = zeros(nx * N, nx);
        Gamma = zeros(nx * N, nu * N);

        A_power = eye(nx);

        for i = 1:N
            A_power = A_power * A;
            rows = (i-1)*nx + (1:nx);
            Phi(rows, :) = A_power;

            for j = 1:i
                cols = (j-1)*nu + (1:nu);
                Gamma(rows, cols) = A^(i-j) * B;
            end
        end

        %% Cost in deviation variables
        Q_b = kron(eye(N), Q);
        R_b = kron(eye(N), R);

        DeltaX_ref_stack = repmat(delta_x_ref, N, 1);

        E = Phi * delta_x0 - DeltaX_ref_stack;

        H = 2 * (Gamma' * Q_b * Gamma + R_b);
        H = (H + H') / 2;
        f = 2 * (Gamma' * Q_b * E);

        %% Rate constraints on delta-U sequence
        D = kron(eye(N), eye(nu)) - kron(diag(ones(N-1,1), -1), eye(nu));

        A_rate = [ D;
                  -D];

        b_rate = [repmat(du_max, N, 1);
                  repmat(du_max, N, 1)];

        %% Absolute input bounds converted to deviation bounds
        % n_min <= u_lin + du_k <= n_max
        lb = repmat([n_min; n_min] - u_lin, N, 1);
        ub = repmat([n_max; n_max] - u_lin, N, 1);

        %% Active circular funnel
        %
        % target_idx = 2 -> pathIds(1)
        % target_idx = 3 -> pathIds(2)
        % ...
        active_funnel_idx = min(max(target_idx - 1, 1), length(pathIds));
        active_node_idx   = pathIds(active_funnel_idx);
        active_poly       = nodes(active_node_idx).poly;

        % Funnel center
        if isfield(nodes, 'c') && ~isempty(nodes(active_node_idx).c)
            funnel_center = nodes(active_node_idx).c(:);
        else
            [cx, cy] = centroid(active_poly);
            funnel_center = [cx; cy];
        end

        % Funnel radius
        % If you have an exact stored radius field, use that instead.
        % Example:
        % funnel_radius = nodes(active_node_idx).R;
        funnel_radius = estimate_circle_radius_from_poly(active_poly, funnel_center);

        %% Initial guess for NLP
        % First try the zero correction sequence, clipped to bounds.
        DU0 = zeros(nu*N, 1);
        DU0 = min(max(DU0, lb), ub);

        % If you want, you can warm-start from previous solution instead.

        %% Objective and nonlinear circle constraints
        objfun = @(DU) 0.5 * DU' * H * DU + f' * DU;

        nonlcon = @(DU) circle_funnel_constraints( ...
            DU, Phi, Gamma, x_lin, delta_x0, funnel_center, funnel_radius, nx, N);

        %% Solve NLP
        [DU_opt, ~, exitflag] = fmincon(objfun, DU0, A_rate, b_rate, ...
            [], [], lb, ub, nonlcon, nlp_options);

        %% Apply input safely
        if exitflag > 0 && ~isempty(DU_opt)
            du_cmd = DU_opt(1:nu);
            u_apply = u_lin + du_cmd;
        else
            warning('fmincon failed at step %d, exitflag=%d. Holding previous input.', t, exitflag);
            u_apply = u_prev;
        end

        u_prev = u_apply;
    end

    %% Plant update
    [xdot, ~] = tugboat3d(x_robot, u_apply);
    x_robot = x_robot + dt_sim * xdot;
    x_robot(3) = atan2(sin(x_robot(3)), cos(x_robot(3)));

    %% Log
    history_x(:, t) = x_robot;
    history_u(:, t) = u_apply;

    %% Stop condition
    if norm(x_robot(1:2) - q_goal') < goal_point_tol
        history_x = history_x(:, 1:t);
        history_u = history_u(:, 1:t);
        break;
    end
end

%% Plotting
% Define figure dimensions to fix "Unrecognized function" error
figWidth = 15; 
figHeight = 12;

% --- Ana Harita ve Rota Grafiği ---
fig1 = figure('WindowState','maximized', 'Color','w'); 
ax1 = gca; hold on; axis equal;
xlim([W(1) W(2)]); ylim([W(3) W(4)]);
set(ax1, 'XTick', [], 'YTick', [], 'Box', 'on', 'LineWidth', 1); 

% Remove Margins for Main Map
set(ax1, 'LooseInset', [0,0,0,0]);
ax1.Position = [0 0 1 1]; 

for i = 1:numel(obs)
    plot(obs{i}, 'FaceColor',[0 0 0], 'FaceAlpha',0.6, 'EdgeColor','none');
end

if ~isempty(pathIds)
    for k = 1:length(pathIds)
        node_idx = pathIds(k);
        plot(nodes(node_idx).poly, 'FaceColor',[1.0 0.85 0.7], 'FaceAlpha',0.40, ...
            'EdgeColor',[1.0 0.5 0.0], 'LineWidth',1.5);
    end
    plot(waypoints(:,1), waypoints(:,2), 'bo', 'MarkerSize',6, 'LineWidth',1.5);
end

plot(history_x(1,:), history_x(2,:), 'r-', 'LineWidth', 2);
plot(q_start(1), q_start(2), 'go', 'MarkerSize',9, 'LineWidth',2);
text(q_start(1), q_start(2), '  start', 'FontSize', 12, 'FontWeight', 'bold');
plot(q_goal(1),  q_goal(2),  'ro', 'MarkerSize',9, 'LineWidth',2);
text(q_goal(1), q_goal(2), '  goal', 'FontSize', 12, 'FontWeight', 'bold');

% --- Kontrol Girişleri Grafiği ---
t_vec = (0:size(history_u,2)-1) * dt_sim;
% Fixed the position logic to use the defined variables
fig2 = figure('Units','centimeters', 'Position',[2, 2, figWidth, figHeight], 'Color','w');
tlo2 = tiledlayout(2,1,'Padding','tight','TileSpacing','compact'); 

nexttile;
plot(t_vec, history_u(1,:), 'b', 'LineWidth',1.5);
grid on; set(gca, 'XTickLabel', [], 'Box', 'on');
ylabel('F_{left} (N)'); 

nexttile;
plot(t_vec, history_u(2,:), 'r', 'LineWidth',1.5);
grid on; set(gca, 'Box', 'on');
xlabel('Time [s]');
ylabel('F_{right} (N)');

% --- Durum Değişkenleri Grafiği ---
fig3 = figure('Units','centimeters', 'Position',[figWidth+3, 2, figWidth, figHeight], 'Color','w');
time_axis = (0:size(history_x,2)-1)*dt_sim;
tlo3 = tiledlayout(3,1,'Padding','tight','TileSpacing','compact');

nexttile;
plot(time_axis, history_x(3,:), 'k', 'LineWidth',1.5);
grid on; set(gca, 'XTickLabel', [], 'YTickLabel', [], 'Box', 'on');
ylabel('Heading');

nexttile;
plot(time_axis, history_x(4,:), 'b', 'LineWidth',1.5); hold on;
plot(time_axis, history_x(5,:), 'r', 'LineWidth',1.5);
grid on; set(gca, 'XTickLabel', [], 'YTickLabel', [], 'Box', 'on');
ylabel('Velocities');

nexttile;
plot(time_axis, history_x(6,:), 'm', 'LineWidth',1.5);
grid on; set(gca, 'XTickLabel', [], 'YTickLabel', [], 'Box', 'on');
ylabel('Yaw Rate');
xlabel('Time [s]');

%% Nonlinear circular funnel constraint
function [c, ceq] = circle_funnel_constraints(DU, Phi, Gamma, x_lin, delta_x0, funnel_center, funnel_radius, nx, N)

    c = zeros(N,1);
    ceq = [];

    for k = 1:N
        rows_k = (k-1)*nx + (1:nx);

        Phi_k   = Phi(rows_k, :);
        Gamma_k = Gamma(rows_k, :);

        xk = x_lin + Phi_k * delta_x0 + Gamma_k * DU;
        pk = xk(1:2);

        c(k) = norm(pk - funnel_center) - funnel_radius;
    end
end

%% Estimate circle radius from stored polygon
function R = estimate_circle_radius_from_poly(polyin, center)

    [vx, vy] = boundary(polyin);

    vx = vx(:);
    vy = vy(:);

    valid = ~(isnan(vx) | isnan(vy));
    vx = vx(valid);
    vy = vy(valid);

    V = [vx vy];

    if isempty(V)
        error('Could not extract vertices from funnel polyshape.');
    end

    if size(V,1) >= 2 && norm(V(1,:) - V(end,:)) < 1e-12
        V(end,:) = [];
    end

    d = sqrt((V(:,1) - center(1)).^2 + (V(:,2) - center(2)).^2);

    % For circular funnels stored as high-edge polygons, these distances
    % should be nearly constant. Mean is a robust estimate.
    R = mean(d);
end