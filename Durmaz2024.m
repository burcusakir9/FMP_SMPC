%% DURMAZ, OZDEMIR & ANKARALI (2024) CONTROLLER + SIMULATION
% "Feedback motion planning via sequential composition of random
%  elliptical funnels" -- circular-funnel special case (a = 1), Eq. (33).
%
% Reduces to the Ege & Ankarali (2019) policy plus one extra
% feedback-linearizing term on omega, (v/rho)*sin(alpha), which is what
% the paper's Proposition 1 uses to guarantee rho is non-increasing
% (i.e. the vehicle provably never leaves the active funnel).

%% CONTROLLER PARAMETERS

Kv       = 0.10;
Ka       = 0.30;
theta0   = 0.0;      % initial heading [rad]

dt_sim   = 0.01;     % simulation sample time [s]
sim_time = 350.0;    % maximum simulation duration [s]
goal_tol = 0.05;     % final stop tolerance on position [m]

% Arrival threshold on rho used only inside the goal funnel, to avoid the
% v/rho singularity at the funnel center (Sec. 3.2.2 "Experimental tuning").
rho_arrival_tol = 0.05;

%% INITIAL STATE

% state = [x; y; theta]
state = [q_start(1); q_start(2); theta0];

time = 0:dt_sim:sim_time;
N = numel(time);

state_hist        = zeros(N, 3);
v_hist            = zeros(N, 1);
omega_hist        = zeros(N, 1);
rho_hist          = zeros(N, 1);
alpha_hist        = zeros(N, 1);
active_funnel_hist = zeros(N, 1);

state_hist(1,:) = state.';

%% SIMULATION

last_idx = N;

for k = 1:N-1

    position = state(1:2).';

    % Select highest-priority funnel containing current vehicle position.
    % pathIds is ordered from start-side funnel toward master/goal funnel.
    [active_path_idx, active_node_id] = selectActiveFunnel(position, nodes, pathIds);

    center = nodes(active_node_id).c;
    is_goal_funnel = (active_path_idx == numel(pathIds));

    % Position relative to active funnel center (local W frame, a = 1)
    x = state(1) - center(1);
    y = state(2) - center(2);

    rho = hypot(x, y);

    % Bearing toward funnel center, Eq. (11) with a = 1
    phi = atan2(-y, -x);

    % Heading error relative to bearing-to-center, Eq. (17)
    alpha = wrapToPiLocal(phi - state(3));

    % Circular-funnel control policy, Eq. (33)
    if is_goal_funnel && rho <= rho_arrival_tol
        % Arrival handling: sidesteps the v/rho singularity at the
        % funnel center once the vehicle is essentially at the goal.
        v = 0;
        omega = 0;
    else
        v     = 2 * Kv * rho * cos(alpha);
        omega = Ka * alpha + (v / rho) * sin(alpha);
    end

    % First-order unicycle model, Eq. (3)
    state_dot = [ ...
        v * cos(state(3));
        v * sin(state(3));
        omega ];

    % Euler integration
    state = state + dt_sim * state_dot;
    state(3) = wrapToPiLocal(state(3));

    % Save
    state_hist(k+1,:)          = state.';
    v_hist(k)                  = v;
    omega_hist(k)               = omega;
    rho_hist(k)                 = rho;
    alpha_hist(k)                = alpha;
    active_funnel_hist(k)      = active_path_idx;

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
rho_hist = rho_hist(1:last_idx);
alpha_hist = alpha_hist(1:last_idx);
active_funnel_hist = active_funnel_hist(1:last_idx);

% Fill final samples for clean plotting
if last_idx > 1
    v_hist(end) = v_hist(end-1);
    omega_hist(end) = omega_hist(end-1);
    rho_hist(end) = rho_hist(end-1);
    alpha_hist(end) = alpha_hist(end-1);
    active_funnel_hist(end) = active_funnel_hist(end-1);
end

fprintf('\n--- Durmaz, Ozdemir & Ankarali (2024) Controller Simulation (circular funnels, a=1) ---\n');
fprintf('Kv           = %.3f\n', Kv);
fprintf('Ka           = %.3f\n', Ka);
fprintf('Simulation   = %.2f s\n', time(end));
fprintf('Final error  = %.4f m\n', norm(state_hist(end,1:2) - q_goal));

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

title(ax1, 'Funnel Chain + Closed-Loop Trajectory (Durmaz et al. 2024, a=1)');

%% PLOT 2: POSITION

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

%% PLOT 3: HEADING

figure('Color','w');

plot(time, rad2deg(state_hist(:,3)), ...
    'LineWidth',1.5);

grid on;

xlabel('Time [s]');
ylabel('\theta [deg]');

title('USV Heading');

%% PLOT 4: POLAR STATES

% figure('Color','w');
%
% yyaxis left
% plot(time, rho_hist, 'LineWidth',1.5);
% ylabel('\rho [m]');
%
% yyaxis right
% plot(time, rad2deg(alpha_hist), 'LineWidth',1.5);
% ylabel('\alpha [deg]');
%
% grid on;
%
% xlabel('Time [s]');
%
% title('Polar States Relative to Active Funnel Center');

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

%% PLOT 6: ACTIVE FUNNEL

% figure('Color','w');
%
% stairs(time, active_funnel_hist, ...
%     'LineWidth',1.5);
%
% grid on;
%
% xlabel('Time [s]');
% ylabel('Active funnel index');
%
% yticks(1:numel(pathIds));
%
% title('Sequential Funnel Switching');

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
