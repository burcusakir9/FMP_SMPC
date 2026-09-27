%% DURMAZ, OZDEMIR & ANKARALI (2024) CONTROLLER + CBF FUNNEL-BOUNDARY FILTER
%
% Durmaz2024.m's nominal circular-funnel navigation law (Eq. 33), UNCHANGED,
% wrapped every step in a minimally-invasive CBF-QP safety filter that
% GUARANTEES the vessel stays inside the currently active funnel circle,
% following Wang, Xiao, Gonzalez-Garcia, Swevers, Ratti & Rus, "Robust
% Model Predictive Control with Control Barrier Functions for Autonomous
% Surface Vessels" (ICRA 2024): they pipe their MPC's nominal control
% through a QP with CBF safety constraints; here the Durmaz2024 law plays
% the role of "nominal control" instead of the MPC.
%
% Barrier (active funnel, center c, radius R):
%   b(x) = R^2 - rho^2,   rho = ||position - c||     (b >= 0 means inside)
%
% For the kinematic unicycle (x=[X;Y;psi], u=[v;omega], f(x)=0, x,y taken
% relative to c):
%
%   b_dot = Lg_b * u,   Lg_b = [-2x*cos(psi) - 2y*sin(psi), 0]
%
% omega does not appear (L_g b w.r.t. omega is zero) -- relative degree 1,
% in v only -- so this is the m=1 ordinary-CBF case of the paper's Def. 3.
% Requiring b_dot + gamma_cbf*b >= 0 gives the QP:
%
%   min ||u - u_nom||^2   s.t.   -Lgb_v * v <= gamma_cbf * b
%
% solved below each step via quadprog. Only v is ever corrected; omega
% (turning) is left entirely to the nominal law. Since Durmaz2024's own
% law already drives rho down under nominal (undisturbed) conditions, the
% filter is expected to sit idle in normal operation -- it is a guarantee
% against violations (e.g. from disturbances or model mismatch), not a
% behavior change by itself.

%% CONTROLLER PARAMETERS (Durmaz2024.m's own nominal law, unchanged)

Kv       = 0.05;
Ka       = 0.30;
theta0   = 0.0;      % initial heading [rad]

dt_sim   = 0.01;     % simulation sample time [s]
sim_time = 350.0;    % maximum simulation duration [s]
goal_tol = 0.05;     % final stop tolerance on position [m]

% Arrival threshold on rho used only inside the goal funnel, to avoid the v/rho singularity at the funnel center
rho_arrival_tol = 0.05;

%% CBF SAFETY-FILTER PARAMETERS

gamma_cbf   = 1.0;   % class-K gain, alpha(b) = gamma_cbf * b
cbf_enabled = true;  % false => bypass the QP, pure nominal Durmaz2024 law

qp_options = optimoptions('quadprog', 'Display', 'off');

%% INITIAL STATE

% state = [x; y; theta]
state = [q_start(1); q_start(2); theta0];

time = 0:dt_sim:sim_time;
N = numel(time);

state_hist         = zeros(N, 3);
v_hist             = zeros(N, 1);
omega_hist         = zeros(N, 1);
vnom_hist          = zeros(N, 1); % nominal (pre-QP) surge command, for comparison
rho_hist           = zeros(N, 1);
alpha_hist         = zeros(N, 1);
active_funnel_hist = zeros(N, 1);
b_hist             = zeros(N, 1); % barrier value R^2 - rho^2 (>=0 means inside the funnel)
funnel_violation   = false(N, 1);
overshoot_hist     = zeros(N, 1); % how far rho exceeds the active funnel radius, if at all
cbf_active_hist    = false(N, 1); % true whenever the QP actually altered v

state_hist(1,:) = state.';

%% SIMULATION

last_idx = N;

for k = 1:N-1

    position = state(1:2).';
    psi = state(3);

    % Select highest-priority funnel containing current vehicle position.
    % pathIds is ordered from start-side funnel toward master/goal funnel.
    [active_path_idx, active_node_id] = selectActiveFunnel(position, nodes, pathIds);

    center = nodes(active_node_id).c;
    R_active = nodes(active_node_id).radius;
    is_goal_funnel = (active_path_idx == numel(pathIds));

    % Position relative to active funnel center
    x = state(1) - center(1);
    y = state(2) - center(2);

    rho = hypot(x, y);

    funnel_violation(k) = rho > R_active;
    overshoot_hist(k) = max(0, rho - R_active);
    b = R_active^2 - rho^2;
    b_hist(k) = b;

    % Bearing toward funnel center
    phi = atan2(-y, -x);

    % Heading error relative to bearing-to-center
    alpha = wrapToPiLocal(phi - psi);

    % Circular-funnel nominal control policy
    if is_goal_funnel && rho <= rho_arrival_tol
        v_nom = 0;
        omega_nom = 0;
    else
        v_nom     = 2 * Kv * rho * cos(alpha);
        omega_nom = Ka * alpha + (v_nom / rho) * sin(alpha);
    end

    % CBF-QP safety filter: min ||u - u_nom||^2 s.t. -Lgb_v*v <= gamma*b.
    if cbf_enabled
        Lgb_v = -2*x*cos(psi) - 2*y*sin(psi);

        H = 2*eye(2);
        f = -2*[v_nom; omega_nom];
        A = [-Lgb_v, 0];
        bineq = gamma_cbf*b;

        [u_opt, ~, exitflag] = quadprog(H, f, A, bineq, [], [], [], [], [v_nom; omega_nom], qp_options);

        if exitflag > 0
            v = u_opt(1);
            omega = u_opt(2);
        else
            % QP infeasible (can happen if already outside the funnel and
            % Lgb_v == 0): fall back to the nominal command rather than
            % fail silently.
            v = v_nom;
            omega = omega_nom;
        end
    else
        v = v_nom;
        omega = omega_nom;
    end

    % 1e-4 m/s, not a tighter tolerance: quadprog's own solver tolerance
    % already perturbs the unconstrained solution by ~1e-8, which would
    % otherwise register as "active" on every step even when the
    % constraint never binds.
    cbf_active_hist(k) = abs(v - v_nom) > 1e-4;

    % First-order unicycle model, Eq. (3)
    state_dot = [ ...
        v * cos(psi);
        v * sin(psi);
        omega ];

    % Euler integration
    state = state + dt_sim * state_dot;
    state(3) = wrapToPiLocal(state(3));

    % Save
    state_hist(k+1,:)          = state.';
    v_hist(k)                  = v;
    omega_hist(k)               = omega;
    vnom_hist(k)                 = v_nom;
    rho_hist(k)                   = rho;
    alpha_hist(k)                  = alpha;
    active_funnel_hist(k)        = active_path_idx;

    % Stop when final goal is reached
    if norm(state(1:2).' - q_goal) <= goal_tol
        last_idx = k + 1;
        break;
    end
end

%% TRIM LOGS

time = time(1:last_idx);
state_hist = state_hist(1:last_idx,:);

v_hist = v_hist(1:last_idx);
omega_hist = omega_hist(1:last_idx);
vnom_hist = vnom_hist(1:last_idx);
rho_hist = rho_hist(1:last_idx);
alpha_hist = alpha_hist(1:last_idx);
active_funnel_hist = active_funnel_hist(1:last_idx);
b_hist = b_hist(1:last_idx);
funnel_violation = funnel_violation(1:last_idx);
overshoot_hist = overshoot_hist(1:last_idx);
cbf_active_hist = cbf_active_hist(1:last_idx);

% Fill final samples for clean plotting
if last_idx > 1
    v_hist(end) = v_hist(end-1);
    omega_hist(end) = omega_hist(end-1);
    vnom_hist(end) = vnom_hist(end-1);
    rho_hist(end) = rho_hist(end-1);
    alpha_hist(end) = alpha_hist(end-1);
    active_funnel_hist(end) = active_funnel_hist(end-1);
    b_hist(end) = b_hist(end-1);
    cbf_active_hist(end) = cbf_active_hist(end-1);
end

fprintf('\n--- Durmaz, Ozdemir & Ankarali (2024) Controller + CBF Funnel Filter ---\n');
fprintf('Kv           = %.3f\n', Kv);
fprintf('Ka           = %.3f\n', Ka);
fprintf('gamma_cbf    = %.3f\n', gamma_cbf);
fprintf('Simulation   = %.2f s\n', time(end));
fprintf('Final error  = %.4f m\n', norm(state_hist(end,1:2) - q_goal));
fprintf('Funnel boundary violations: %d / %d steps (%.2f%%)\n', ...
    sum(funnel_violation), numel(funnel_violation), 100*mean(funnel_violation));
fprintf('Max overshoot beyond active funnel radius: %.4f m\n', max(overshoot_hist));
fprintf('CBF actively correcting v: %d / %d steps (%.2f%%)\n', ...
    sum(cbf_active_hist), numel(cbf_active_hist), 100*mean(cbf_active_hist));

%% PLOT 1: CLOSED-LOOP TRAJECTORY

fig1 = figure('WindowState','maximized', 'Color','w');
ax1 = axes('Parent', fig1);
hold(ax1, 'on');
axis(ax1, 'equal');

xlim(ax1, [W(1) W(2)]);
ylim(ax1, [W(3) W(4)]);

set(ax1, 'XTick', [], 'YTick', [], 'Box', 'on');
set(ax1, 'LooseInset', [0,0,0,0]);
ax1.Position = [0 0 1 1];

% Obstacles
for i = 1:numel(obs)
    plot(ax1, obs{i}, ...
        'FaceColor',[0 0 0], ...
        'FaceAlpha',0.6, ...
        'EdgeColor','none');
end

% Funnel chain
for k = 1:numel(pathIds)
    node_id = pathIds(k);

    plot(ax1, nodes(node_id).poly, ...
        'FaceColor',[1.0 0.85 0.7], ...
        'FaceAlpha',0.30, ...
        'EdgeColor',[1.0 0.5 0.0], ...
        'LineWidth',1.2);

    plot(ax1, nodes(node_id).c(1), nodes(node_id).c(2), ...
        '.', 'Color',[0.85 0.35 0.0], 'MarkerSize',12);
end

% Executed trajectory
plot(ax1, state_hist(:,1), state_hist(:,2), ...
    'b-', 'LineWidth',2.0);

% Start and goal
plot(ax1, q_start(1), q_start(2), ...
    'go', 'MarkerSize',9, 'LineWidth',2);

text(ax1, q_start(1), q_start(2), ...
    '  start', 'FontSize',12, 'FontWeight','bold');

plot(ax1, q_goal(1), q_goal(2), ...
    'ro', 'MarkerSize',9, 'LineWidth',2);

text(ax1, q_goal(1), q_goal(2), ...
    '  goal', 'FontSize',12, 'FontWeight','bold');

title(ax1, 'Funnel Chain + Closed-Loop Trajectory (Durmaz2024 + CBF Funnel Filter)');

%% PLOT 2: BARRIER VALUE AND CBF-FILTERED VS NOMINAL SURGE

figure('Color','w');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile;
plot(time, b_hist, 'LineWidth',1.5);
hold on;
yline(0, 'r--', 'LineWidth',1.2);
grid on; ylabel('b = R^2 - \rho^2'); xlabel('Time [s]');
title('Barrier Value (b < 0 means outside the active funnel)');

nexttile;
plot(time, vnom_hist, 'LineWidth',1.2); hold on;
plot(time, v_hist, 'LineWidth',1.2);
grid on; ylabel('v [m/s]'); xlabel('Time [s]');
legend('v_{nominal} (Durmaz2024)','v_{CBF-filtered}','Location','best');
title('Nominal vs Safety-Filtered Surge Command');

%% PLOT 3: POSITION

figure('Color','w');

plot(time, state_hist(:,1), 'LineWidth',1.5);
hold on;
plot(time, state_hist(:,2), 'LineWidth',1.5);

yline(q_goal(1), '--', 'LineWidth',1.0);
yline(q_goal(2), '--', 'LineWidth',1.0);

grid on;

xlabel('Time [s]');
ylabel('Position [m]');

legend('x', 'y', 'x_{goal}', 'y_{goal}', ...
    'Location','best');

title('USV Position');

%% PLOT 4: HEADING

figure('Color','w');

plot(time, rad2deg(state_hist(:,3)), ...
    'LineWidth',1.5);

grid on;

xlabel('Time [s]');
ylabel('\theta [deg]');

title('USV Heading');

%% PLOT 5: CONTROL INPUTS

figure('Color','w');

yyaxis left
plot(time, v_hist, 'LineWidth',1.5);
ylabel('v [m/s]');

yyaxis right
plot(time, omega_hist, 'LineWidth',1.5);
ylabel('\omega [rad/s]');

grid on;

xlabel('Time [s]');

title('Control Inputs');

%% FUNCTIONS

function [active_path_idx, active_node_id] = selectActiveFunnel(position, nodes, pathIds)

    % Highest priority is the funnel closest to the master/goal funnel.
    % Since pathIds = [start-side ... goal-side], search backwards.

    active_path_idx = 1;
    active_node_id = pathIds(1);

    for k = numel(pathIds):-1:1

        node_id = pathIds(k);

        distance_to_center = norm(position - nodes(node_id).c);

        if distance_to_center <= nodes(node_id).radius
            active_path_idx = k;
            active_node_id = node_id;
            return;
        end
    end
end

function angle = wrapToPiLocal(angle)

    angle = mod(angle + pi, 2*pi) - pi;

end
