%% Parameters

dt_sim = 0.02;          % Simulation sample time
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
dn_max = 200;
du_max = [dn_max; dn_max];

% Waypoint / stopping tolerances
wp_tol = 2.0;
goal_point_tol = 0.5;

% Stage cost weights
Qp = diag([40, 40, 15]);       % [x y yaw]
Qdyn = diag([0, 2, 1]);        % [u v r]  -> do NOT penalize surge u
R = diag([0.002, 0.002]);      % lower input penalty
Rd = diag([0.02, 0.02]);       % rate penalty

% Terminal cost (position/yaw terminal weight)
Qf = diag([120, 120, 40]);

% Desired nominal forward speed toward waypoint
u_ref_mag = 1.0;               % can tune: 0.5 to 1.5

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
    target_idx = 2;
else
    target_idx = 1;
end

history_x = zeros(nx, T);
history_u = zeros(nu, T);

u_prev = [0; 0];
u_apply = [0; 0];

fprintf('Starting tugboat nonlinear MPC waypoint tracking.\n');

U_guess = zeros(nu * N, 1); %#ok<NASGU>

%% LTV-MPC Simulation Loop
for t = 1:T
    if mod(t-1, mpc_interval) == 0
        % 1. Waypoint & Reference Logic
        target = waypoints(target_idx, :)';

        if norm(x_robot(1:2) - target) < wp_tol && target_idx < size(waypoints, 1)
            target_idx = target_idx + 1;
            target = waypoints(target_idx, :)';
        end
        
        dvec = target - x_robot(1:2);
        yaw_ref_raw = atan2(dvec(2), dvec(1));
        yaw_err = atan2(sin(yaw_ref_raw - x_robot(3)), cos(yaw_ref_raw - x_robot(3)));
        yaw_ref = x_robot(3) + yaw_err;
        
        dist_to_wp = norm(target - x_robot(1:2));
        u_ref = min(u_ref_mag, 0.5 * dist_to_wp);
        x_ref = [target(1); target(2); yaw_ref; u_ref; 0; 0];

        % Linearize dynamics
        [A_sys, B_sys] = tugboat3d_linearized(x_robot, u_prev);
        
        % Discretize
        Ad = eye(nx) + A_sys * dt_mpc;
        Bd = B_sys * dt_mpc;
        
        % Build QP matrices
        [H, f, A_cons, b_cons] = build_qp_matrices(Ad, Bd, x_robot, x_ref, u_prev, ...
                                     N, Qp, Qdyn, Qf, R, Rd, n_min, n_max, du_max);
        
        % Solve QP for optimal control sequence U (not dU)
        qp_options = optimoptions('quadprog', ...
            'Display', 'off', ...
            'Algorithm', 'interior-point-convex');

        [U_opt, ~, exitflag] = quadprog(H, f, A_cons, b_cons, [], [], [], [], [], qp_options);
        
        if exitflag > 0
            % Apply first control input
            u_apply = U_opt(1:nu);
            u_prev = u_apply;
        else
            % Fallback
            u_apply = u_prev;
        end
    end

    % 5. Plant Update (Keep nonlinear for reality)
    [xdot, ~] = tugboat3d(x_robot, u_apply);
    x_robot = x_robot + dt_sim * xdot;
    x_robot(3) = atan2(sin(x_robot(3)), cos(x_robot(3)));
    
    % Log
    history_x(:, t) = x_robot;
    history_u(:, t) = u_apply;

    % Stop if final goal reached
    if norm(x_robot(1:2) - q_goal(:)) < goal_point_tol && target_idx == size(waypoints,1)
        history_x = history_x(:, 1:t);
        history_u = history_u(:, 1:t);
        fprintf('Goal reached at t = %.2f s\n', (t-1)*dt_sim);
        break;
    end
end

%% Plotting
figure('WindowState','maximized', 'Color','w'); hold on; axis equal;
xlim([W(1) W(2)]); ylim([W(3) W(4)]);
title('Figür 2: Path, Waypoints, and Tugboat Tracking');

for i = 1:numel(obs)
    plot(obs{i}, 'FaceColor',[0 0 0], 'FaceAlpha',0.6, 'EdgeColor','none');
end

if ~isempty(pathIds)
    for k = 1:length(pathIds)
        node_idx = pathIds(k);
        plot(nodes(node_idx).poly, ...
            'FaceColor',[1.0 0.85 0.7], ...
            'FaceAlpha',0.40, ...
            'EdgeColor',[1.0 0.5 0.0], ...
            'LineWidth',1.5);
    end

    plot(waypoints(:,1), waypoints(:,2), 'bo', 'MarkerSize',6, 'LineWidth',1.5);
end

plot(q_start(1), q_start(2), 'go', 'MarkerSize',9, 'LineWidth',2);
plot(q_goal(1),  q_goal(2),  'ro', 'MarkerSize',9, 'LineWidth',2);
plot(history_x(1,:), history_x(2,:), 'r-', 'LineWidth', 2);

grid on;
xlabel('X [m]');
ylabel('Y [m]');
legend('Obstacles / Funnels', 'Waypoints', 'Start', 'Goal', 'Robot Path');

t_vec = (0:size(history_u,2)-1) * dt_sim;

figure('Color','w');
subplot(2,1,1);
plot(t_vec, history_u(1,:), 'b', 'LineWidth',1.5);
grid on;
ylabel('n_1');
title('Propeller Inputs');

subplot(2,1,2);
plot(t_vec, history_u(2,:), 'r', 'LineWidth',1.5);
grid on;
xlabel('Time [s]');
ylabel('n_2');

figure('Color','w');
subplot(3,1,1);
plot((0:size(history_x,2)-1)*dt_sim, history_x(3,:), 'k', 'LineWidth',1.5);
grid on;
ylabel('\psi [rad]');
title('Yaw and Body Velocities');

subplot(3,1,2);
plot((0:size(history_x,2)-1)*dt_sim, history_x(4,:), 'b', 'LineWidth',1.5); hold on;
plot((0:size(history_x,2)-1)*dt_sim, history_x(5,:), 'r', 'LineWidth',1.5);
grid on;
ylabel('u, v [m/s]');
legend('u','v');

subplot(3,1,3);
plot((0:size(history_x,2)-1)*dt_sim, history_x(6,:), 'm', 'LineWidth',1.5);
grid on;
xlabel('Time [s]');
ylabel('r [rad/s]');

%% ===== Helper functions =====

function J = nmpc_cost(U, x0, x_ref, u_prev, N, dt_mpc, Qp, Qdyn, Qf, R, Rd)
    nx = length(x0);
    nu = 2;

    X = zeros(nx, N+1);
    X(:,1) = x0;

    J = 0;
    u_last = u_prev;

    for k = 1:N
        uk = U((k-1)*nu + (1:nu));

        [xdot, ~] = tugboat3d(X(:,k), uk);
        X(:,k+1) = X(:,k) + dt_mpc * xdot;
        X(3,k+1) = atan2(sin(X(3,k+1)), cos(X(3,k+1)));

        ep = [X(1,k+1)-x_ref(1);
              X(2,k+1)-x_ref(2);
              atan2(sin(X(3,k+1)-x_ref(3)), cos(X(3,k+1)-x_ref(3)))];

        ed = [X(4,k+1)-x_ref(4);
              X(5,k+1)-x_ref(5);
              X(6,k+1)-x_ref(6)];

        du = uk - u_last;

        J = J + ep' * Qp * ep + ed' * Qdyn * ed + uk' * R * uk + du' * Rd * du;

        u_last = uk;
    end

    epN = [X(1,N+1)-x_ref(1);
           X(2,N+1)-x_ref(2);
           atan2(sin(X(3,N+1)-x_ref(3)), cos(X(3,N+1)-x_ref(3)))];

    J = J + epN' * Qf * epN;
end

function [c, ceq] = nmpc_constraints(U, x0, u_prev, N, dt_mpc, du_max, funnel_center, funnel_radius)
    nx = length(x0);
    nu = 2;

    X = zeros(nx, N+1);
    X(:,1) = x0;

    c_rate = zeros(2*nu*N, 1);
    c_funnel = zeros(N, 1);

    u_last = u_prev;

    for k = 1:N
        uk = U((k-1)*nu + (1:nu));

        du = uk - u_last;
        idx = (k-1)*2*nu + (1:2*nu);

        c_rate(idx(1:nu)) = du - du_max;
        c_rate(idx(nu+1:end)) = -du - du_max;

        [xdot, ~] = tugboat3d(X(:,k), uk);
        X(:,k+1) = X(:,k) + dt_mpc * xdot;
        X(3,k+1) = atan2(sin(X(3,k+1)), cos(X(3,k+1)));

        px = X(1,k+1);
        py = X(2,k+1);
        dx = px - funnel_center(1);
        dy = py - funnel_center(2);

        c_funnel(k) = dx^2 + dy^2 - funnel_radius^2;

        u_last = uk;
    end

    c = [c_rate;
         c_funnel];
    ceq = [];
end

function R_est = estimateNodeRadius(node)
    if isfield(node, 'r') && ~isempty(node.r)
        R_est = node.r;
        return;
    end

    c = node.c(:)';
    [vx, vy] = boundary(node.poly);

    if isempty(vx)
        R_est = 1.0;
        return;
    end

    d = sqrt((vx - c(1)).^2 + (vy - c(2)).^2);
    R_est = min(d);

    if isempty(R_est) || R_est <= 0
        R_est = 1.0;
    end
end

function [H, f, A_ineq, b_ineq] = build_qp_matrices(Ad, Bd, x0, x_ref, u_prev, N, Qp, Qdyn, Qf, R, Rd, n_min, n_max, du_max)
    nx = 6;
    nu = 2;

    % Stage and terminal state weights
    Q_block = blkdiag(Qp, Qdyn);      % 6x6
    Qf_full = blkdiag(Qf, zeros(3));  % 6x6

    % Prediction matrices:
    % X = Sx*x0 + Su*U
    % X = [x1; x2; ...; xN], U = [u0; u1; ...; u_{N-1}]
    Sx = zeros(nx*N, nx);
    Su = zeros(nx*N, nu*N);

    for i = 1:N
        Sx((i-1)*nx+1:i*nx, :) = Ad^i;

        for j = 1:i
            Su((i-1)*nx+1:i*nx, (j-1)*nu+1:j*nu) = Ad^(i-j) * Bd;
        end
    end

    % State cost matrix
    Q_cells = cell(N,1);
    for i = 1:N-1
        Q_cells{i} = Q_block;
    end
    Q_cells{N} = Qf_full;
    Q_bar = blkdiag(Q_cells{:});      % 6N x 6N

    % Input effort cost
    R_bar = kron(eye(N), R);          % 2N x 2N

    % Reference stack
    X_ref = repmat(x_ref, N, 1);      % 6N x 1

    % Free response error
    E_free = Sx * x0 - X_ref;         % 6N x 1

    % Difference matrix for input increments:
    % dU = D*U - d0
    % dU(1) = u0 - u_prev
    % dU(k) = uk - u_{k-1}, k>=2
    D = zeros(nu*N, nu*N);
    for k = 1:N
        rows = (k-1)*nu+1:k*nu;
        cols = (k-1)*nu+1:k*nu;
        D(rows, cols) = eye(nu);

        if k >= 2
            prev_cols = (k-2)*nu+1:(k-1)*nu;
            D(rows, prev_cols) = -eye(nu);
        end
    end

    d0 = [u_prev; zeros(nu*(N-1),1)];

    % Rate cost
    Rd_bar = kron(eye(N), Rd);

    % Hessian
    H = 2 * (Su' * Q_bar * Su + R_bar + D' * Rd_bar * D);
    H = (H + H') / 2;
    H = H + 1e-8 * eye(size(H));

    % Gradient
    f = 2 * (Su' * Q_bar * E_free - D' * Rd_bar * d0);

    % ---------- Inequality constraints ----------

    % Make sure bounds are 2x1 vectors
    if isscalar(n_max)
        n_max = [n_max; n_max];
    end
    if isscalar(n_min)
        n_min = [n_min; n_min];
    end
    if isscalar(du_max)
        du_max = [du_max; du_max];
    end

    % 1) Rate constraints:
    % -du_max <= D*U - d0 <= du_max
    du_stack = repmat(du_max(:), N, 1);   % 2N x 1

    A_rate = [ D;
              -D];                        % 4N x 2N
    b_rate = [ du_stack + d0;
               du_stack - d0];            % 4N x 1

    % 2) Input bounds:
    % n_min <= uk <= n_max  for all k
    u_max_stack = repmat(n_max(:), N, 1); % 2N x 1
    u_min_stack = repmat(n_min(:), N, 1); % 2N x 1

    A_input = [ eye(nu*N);
               -eye(nu*N)];               % 4N x 2N
    b_input = [ u_max_stack;
               -u_min_stack];             % 4N x 1

    % Final inequality system
    A_ineq = [A_rate;
              A_input];                   % 8N x 2N
    b_ineq = [b_rate;
              b_input];                   % 8N x 1

    % Cleanup
    A_ineq(isnan(A_ineq)) = 0;
    A_ineq(isinf(A_ineq)) = 0;
    b_ineq(isnan(b_ineq)) = 0;
    b_ineq(isinf(b_ineq)) = 1e12;
end