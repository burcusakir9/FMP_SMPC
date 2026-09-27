%% RSC CIRCULAR FUNNELS + SWAY-AWARE PARTIAL FEEDBACK LINEARIZATION
%
% Proposed extension of the Ege2019 / Durmaz2024 circular-funnel control
% laws from a kinematic unicycle to the actual 3-DOF dynamic tugboat3d.m
% plant (state = [X;Y;psi;u;v;r;F_L;F_R], first-order actuator lag).
%
% Architecture: RSC funnel -> outer sway-aware (rho,alpha) law -> virtual
% (u_c, r_c) -> inner exact partial feedback linearization -> (X, N) ->
% thrust allocation -> (F_L_cmd, F_R_cmd).
%
% Outer loop (exact local kinematics with sway, eta_dot = R(psi)*nu):
%   rho    = ||position - funnel_center||
%   alpha  = bearing_to_center - psi
%   rho_dot   = -u*cos(alpha) - v*sin(alpha)
%   alpha_dot = (u*sin(alpha) - v*cos(alpha))/rho - r
%
%   u_c = k_rho*rho*cos(alpha)
%   r_c = k_alpha*alpha + (u_c*sin(alpha) - v*cos(alpha))/rho
%
% Substituting u=u_c, r=r_c gives alpha_dot = -k_alpha*alpha exactly (the
% sway-aware analogue of Durmaz2024's feedback-linearizing term).
%
% Inner loop: u_c and r_c are differentiated EXACTLY (true backstepping,
% no dynamic-surface-control command filters and hence no boundary-layer
% error term). u_c_dot is immediate from the exact kinematics above and
% needs no v_dot:
%
%   u_c_dot = k_rho*(rho_dot*cos(alpha) - rho*sin(alpha)*alpha_dot)
%
% r_c_dot is harder: it contains v_dot, and v_dot depends on r_dot, which
% is precisely what N is being solved for. That is a genuine algebraic
% loop -- but a LINEAR one, so it is closed in place rather than dodged.
% Writing g = u_c*sin(alpha) - v*cos(alpha) and r_c = k_alpha*alpha +
% g/rho, the v_dot dependence is affine:
%
%   r_c_dot = A + B*v_dot,        B = -cos(alpha)/rho
%
% Substituting row 2 of the plant, v_dot = (-f2 - m23*r_dot)/m22, turns
% this into r_c_dot = C + D*r_dot with
%
%   C = A - B*f2/m22,   D = -B*m23/m22 = (m23/m22)*cos(alpha)/rho
%
% and imposing the desired error dynamics r_dot = r_c_dot - k_r*(r-r_c)
% gives the closed-form solve
%
%   r_dot_des = (C - k_r*(r - r_c)) / (1 - D)
%
% The only singularity is 1 - D = 0, i.e. cos(alpha)/rho = m22/m23 =
% 5.838, i.e. rho = 0.171 m at worst. rho_floor (below) already bounds
% rho from beneath at 1.0 m, so |D| <= 0.171 and (1 - D) >= 0.83 -- a
% factor of ~5 of margin that never closes. D_max is logged at the end so
% this can be confirmed empirically rather than assumed.
%
% Only row 1 of the vessel's M-matrix is decoupled from sway (X enters
% only the surge equation), so X is an EXACT feedback-linearizing input.
% Rows 2-3 (sway/yaw) are coupled through the off-diagonal added-mass
% term m23, and N alone cannot independently set both v_dot and r_dot --
% v is a genuinely unactuated internal state. Solving row 2 for v_dot and
% substituting into row 3 gives the exact partial-feedback-linearizing N
% (v_dot itself is left to evolve according to the true coupled dynamics,
% i.e. collocated partial linearization, not the diagonal-model
% approximation):
%
%   X = m11*(u_c_dot - k_u*(u-u_c)) + f1(nu)
%   N = (m33 - m23^2/m22)*r_dot_des - (m23/m22)*f2(nu) + f3(nu)
%
% where f1,f2,f3 are the row sums of C(nu)*nu + D(nu)*nu from
% tugboat3d.m, reproduced exactly below so the cancellation matches the
% actual simulated plant.
%
% A close-range instability showed up in testing that is directly
% relevant to the "v(t) -> 0" assumption needed for the asymptotic
% convergence argument: near rho -> 0, sway feeding the v*cos(alpha)/rho
% term in r_c excites more sway through the m23 sway-yaw coupling faster
% than the outer-loop gains damp it, producing a sustained limit cycle
% orbiting the goal instead of settling (this does not occur in the
% kinematic Ege2019/Durmaz2024 controllers, which have no such coupling).
% See the rho_floor comment below for the fix used here.
%
% NOT implemented in this first version (left as documented future
% work, matching steps 8-9 of the proposed development order):
%   - backstepping through the thruster's own first-order lag
%     (F_cmd is commanded equal to the desired thrust directly);
%   - shrinking funnel radii by a velocity-dependent braking margin.
% Instead, funnel-boundary violations (rho exceeding the active funnel's
% radius -- the concern raised about positional-only funnel invariance
% not holding for a vehicle with real momentum) are logged as a
% diagnostic so the effect can be measured empirically.

%% GAINS
%
% Each gain falls back to its default only if not already set in the
% workspace, so a driver script can inject a different gain set (e.g. a
% stress test) without editing this file.

if ~exist('k_rho','var'),   k_rho   = 0.10; end % outer surge command gain [1/s]
if ~exist('k_alpha','var'), k_alpha = 0.30;  end % outer yaw-rate command gain [1/s]

if ~exist('k_u','var'), k_u = 1.0; end % inner surge tracking rate [1/s]
if ~exist('k_r','var'), k_r = 1.0; end % inner yaw-rate tracking rate [1/s]

if ~exist('F_max','var'), F_max = 100;  end % thruster force saturation [N], matches tugboat_mpc.m
if ~exist('F_min','var'), F_min = -100; end

if ~exist('dt_sim','var'),   dt_sim   = 0.01;  end
if ~exist('sim_time','var'), sim_time = 500.0; end
if ~exist('goal_tol','var'), goal_tol = 0.05;  end

if ~exist('rho_floor','var'), rho_floor = 1.0; end
% Floors the (.)/rho term in r_c well before rho reaches 0 -- see the
% close-range instability note above for why this is needed.

%% INITIAL STATE

% state = [X; Y; psi; u; v; r; F_L; F_R]
state = [q_start(1); q_start(2); 0; 0; 0; 0; 0; 0];

time = 0:dt_sim:sim_time;
N_steps = numel(time);

state_hist          = zeros(N_steps, 8);
v_hist               = zeros(N_steps, 1); % vessel speed through water, hypot(u,v)
omega_hist           = zeros(N_steps, 1); % yaw rate r
rho_hist             = zeros(N_steps, 1);
alpha_hist           = zeros(N_steps, 1);
X_hist               = zeros(N_steps, 1);
Ncmd_hist            = zeros(N_steps, 1);
active_funnel_hist   = zeros(N_steps, 1);
funnel_violation     = false(N_steps, 1);
overshoot_hist       = zeros(N_steps, 1); % rho - active_radius, clipped at 0 below
Fcmd_hist            = zeros(N_steps, 2); % UNSATURATED thrust command [F_L F_R]
sat_hist             = false(N_steps, 1); % true when either channel clipped
eu_hist              = zeros(N_steps, 1); % surge command-tracking error u_c - u
er_hist              = zeros(N_steps, 1); % yaw   command-tracking error r_c - r

state_hist(1,:) = state.';

%% SIMULATION

last_idx = N_steps;

rc_D_max = 0; % largest |rc_D| seen; margin to the (1 - D) = 0 singularity

for k = 1:N_steps-1

    position = state(1:2).';
    psi = state(3);
    u   = state(4);
    v   = state(5);
    r   = state(6);

    % Select highest-priority funnel containing current vehicle position.
    [active_path_idx, active_node_id] = selectActiveFunnel(position, nodes, pathIds);

    center = nodes(active_node_id).c;
    R_active = nodes(active_node_id).radius;

    % Relative position to active funnel center
    ex = state(1) - center(1);
    ey = state(2) - center(2);
    rho = hypot(ex, ey);

    % Bearing toward funnel center and sway-aware heading error
    phi   = atan2(-ey, -ex);
    alpha = wrapToPiLocal(phi - psi);

    funnel_violation(k) = rho > R_active;
    overshoot_hist(k) = max(0, rho - R_active);

    % Floored range used everywhere the (.)/rho division appears. When the
    % floor is active rho_e is constant, so its derivative is zero -- that
    % has to be respected below or the differentiation is inconsistent
    % with the expression actually being used.
    rho_e = max(rho, rho_floor);

    % Exact kinematics with sway (eta_dot = R(psi)*nu)
    rho_dot   = -u * cos(alpha) - v * sin(alpha);
    alpha_dot = (u * sin(alpha) - v * cos(alpha)) / rho_e - r;

    rho_e_dot = (rho > rho_floor) * rho_dot;

    % Outer sway-aware funnel control policy. u_c already tapers to 0 as
    % rho -> 0 on its own; only the (.)/rho term in r_c needs the floor.
    uc = k_rho * rho * cos(alpha);
    g  = uc * sin(alpha) - v * cos(alpha);
    rc = k_alpha * alpha + g / rho_e;

    % Exact derivative of the surge command -- no v_dot appears, so this
    % one closes immediately.
    uc_dot = k_rho * (rho_dot * cos(alpha) - rho * sin(alpha) * alpha_dot);

    % Inner-loop tracking errors against the TRUE commands (no filters)
    eu = u - uc;
    er = r - rc;

    % Exact partial feedback linearization using the true plant model
    [f1, f2, f3, m11, m22, m23, m33] = tugboatModelTerms([u; v; r]);

    % Exact derivative of the yaw-rate command. g_dot splits into a part
    % known from measured state and a part proportional to v_dot:
    %   g_dot = g_dot_known - v_dot*cos(alpha)
    g_dot_known = uc_dot * sin(alpha) ...
                + (uc * cos(alpha) + v * sin(alpha)) * alpha_dot;

    % r_c_dot = A + B*v_dot
    rc_A = k_alpha * alpha_dot ...
      + (g_dot_known * rho_e - g * rho_e_dot) / rho_e^2;
    rc_B = -cos(alpha) / rho_e;

    % Close the algebraic loop: v_dot = (-f2 - m23*r_dot)/m22 turns the
    % above into r_c_dot = C + D*r_dot.
    rc_C = rc_A - rc_B * f2 / m22;
    rc_D = -rc_B * m23 / m22;

    rc_D_max = max(rc_D_max, abs(rc_D));

    % Impose r_dot = r_c_dot - k_r*er and solve for r_dot. The guard is
    % defensive only: with rho_floor = 1.0, |D| <= 0.171 always.
    denom = 1 - rc_D;
    if abs(denom) < 0.1
        denom = 0.1 * sign(denom + (denom == 0));
    end

    r_dot_des = (rc_C - k_r * er) / denom;

    X  = m11 * (uc_dot - k_u * eu) + f1;
    Nc = (m33 - m23^2 / m22) * r_dot_des - (m23 / m22) * f2 + f3;

    % Thrust allocation, matches tugboat3d.m's B_prop = [1 1; 0 0; d -d]
    d_beam = 0.29 / 2;
    FL_unsat = (X + Nc / d_beam) / 2;
    FR_unsat = (X - Nc / d_beam) / 2;

    FL_cmd = min(max(FL_unsat, F_min), F_max);
    FR_cmd = min(max(FR_unsat, F_min), F_max);

    % Log the UNSATURATED command -- shows how much control authority
    % the law is actually asking for, not just what the plant received.
    Fcmd_hist(k,:) = [FL_unsat, FR_unsat];
    sat_hist(k)    = (FL_cmd ~= FL_unsat) || (FR_cmd ~= FR_unsat);
    eu_hist(k)     = uc - u;
    er_hist(k)     = rc - r;

    % Plant update: full 8-state dynamic tugboat model, incl. actuator lag
    xdot = tugboat3d(state, [FL_cmd; FR_cmd]);
    state = state + dt_sim * xdot;
    state(3) = wrapToPiLocal(state(3));

    % Save
    state_hist(k+1,:)        = state.';
    v_hist(k)                = hypot(state(4), state(5));
    omega_hist(k)             = state(6);
    rho_hist(k)                = rho;
    alpha_hist(k)                = alpha;
    X_hist(k)                      = X;
    Ncmd_hist(k)                     = Nc;
    active_funnel_hist(k)              = active_path_idx;

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
X_hist = X_hist(1:last_idx);
Ncmd_hist = Ncmd_hist(1:last_idx);
active_funnel_hist = active_funnel_hist(1:last_idx);
funnel_violation = funnel_violation(1:last_idx);
overshoot_hist = overshoot_hist(1:last_idx);
Fcmd_hist = Fcmd_hist(1:last_idx,:);
sat_hist = sat_hist(1:last_idx);
eu_hist = eu_hist(1:last_idx);
er_hist = er_hist(1:last_idx);

if last_idx > 1
    v_hist(end) = v_hist(end-1);
    omega_hist(end) = omega_hist(end-1);
    rho_hist(end) = rho_hist(end-1);
    alpha_hist(end) = alpha_hist(end-1);
    X_hist(end) = X_hist(end-1);
    Ncmd_hist(end) = Ncmd_hist(end-1);
    active_funnel_hist(end) = active_funnel_hist(end-1);
end

fprintf('\n--- RSC + Sway-Aware Partial Feedback Linearization Simulation ---\n');
fprintf('k_rho = %.3f, k_alpha = %.3f, k_u = %.3f, k_r = %.3f\n', k_rho, k_alpha, k_u, k_r);
fprintf('Simulation   = %.2f s\n', time(end));
fprintf('Final error  = %.4f m\n', norm(state_hist(end,1:2) - q_goal));
fprintf('Funnel boundary violations: %d / %d steps (%.2f%%)\n', ...
    sum(funnel_violation), numel(funnel_violation), 100*mean(funnel_violation));
fprintf('Max overshoot beyond active funnel radius: %.4f m\n', max(overshoot_hist));
fprintf('Peak |thrust command|: %.2f N  (limit %.0f N)\n', max(abs(Fcmd_hist(:))), F_max);
fprintf('Steps saturated: %d / %d (%.2f%%)\n', sum(sat_hist), numel(sat_hist), 100*mean(sat_hist));
fprintf('RMS surge tracking error |u_c - u|: %.4f m/s\n', sqrt(mean(eu_hist.^2)));
fprintf('RMS yaw   tracking error |r_c - r|: %.4f rad/s\n', sqrt(mean(er_hist.^2)));
fprintf('Max |D| in algebraic loop: %.4f  (singularity at |D| = 1)\n', rc_D_max);

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

for i = 1:numel(obs)
    plot(ax1, obs{i}, 'FaceColor',[0 0 0], 'FaceAlpha',0.6, 'EdgeColor','none');
end

for k = 1:numel(pathIds)
    node_id = pathIds(k);

    plot(ax1, nodes(node_id).poly, 'FaceColor',[1.0 0.85 0.7], 'FaceAlpha',0.30, ...
        'EdgeColor',[1.0 0.5 0.0], 'LineWidth',1.2);

    plot(ax1, nodes(node_id).c(1), nodes(node_id).c(2), ...
        '.', 'Color',[0.85 0.35 0.0], 'MarkerSize',12);
end

plot(ax1, state_hist(:,1), state_hist(:,2), 'b-', 'LineWidth',2.0);

plot(ax1, q_start(1), q_start(2), 'go', 'MarkerSize',9, 'LineWidth',2);
text(ax1, q_start(1), q_start(2), '  start', 'FontSize',12, 'FontWeight','bold');

plot(ax1, q_goal(1), q_goal(2), 'ro', 'MarkerSize',9, 'LineWidth',2);
text(ax1, q_goal(1), q_goal(2), '  goal', 'FontSize',12, 'FontWeight','bold');

title(ax1, 'Funnel Chain + Closed-Loop Trajectory (Sway-Aware PFL)');

%% PLOT 2: VELOCITIES

figure('Color','w');
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

nexttile;
plot(time, state_hist(:,4), 'LineWidth',1.5); hold on;
plot(time, state_hist(:,5), 'LineWidth',1.5);
grid on; ylabel('[m/s]'); legend('u','v','Location','best');
title('Surge / Sway Velocity');

nexttile;
plot(time, state_hist(:,6), 'LineWidth',1.5);
grid on; ylabel('r [rad/s]');
title('Yaw Rate');

nexttile;
plot(time, rho_hist, 'LineWidth',1.5); hold on;
plot(time, rad2deg(alpha_hist), 'LineWidth',1.5);
grid on; xlabel('Time [s]'); legend('\rho [m]', '\alpha [deg]', 'Location','best');
title('Polar States Relative to Active Funnel Center');

%% PLOT 3: THRUST COMMANDS

figure('Color','w');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile;
plot(time, X_hist, 'LineWidth',1.5); hold on;
plot(time, Ncmd_hist, 'LineWidth',1.5);
grid on; ylabel('Force / Moment'); legend('X [N]','N [N\cdotm]','Location','best');
title('Feedback-Linearizing Commands');

nexttile;
plot(time, state_hist(:,7), 'LineWidth',1.5); hold on;
plot(time, state_hist(:,8), 'LineWidth',1.5);
grid on; ylabel('F [N]'); xlabel('Time [s]'); legend('F_L','F_R','Location','best');
title('Actual Thruster Forces');

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

function [f1, f2, f3, m11, m22, m23, m33] = tugboatModelTerms(nu)

    % Reproduces the exact M, C(nu), D(nu) construction of tugboat3d.m so
    % the dynamic-inversion terms above exactly cancel the true plant.

    m  = 10.2;
    Iz = 0.63994;

    Xu    = -5.76909;
    Xuu   = -2.17161;
    Xudot = -0.87818;

    Yv    = -3.98659;
    Yr    = -0.0001;
    Yvv   = -3.95131;
    Yvdot = -1.05279;
    Yrdot = -1.92760;

    Nv    = -0.0001;
    Nr    = -0.12392;
    Nrr   = -0.33077;
    Nrdot = -0.04531;

    m11 = m - Xudot;
    m22 = m - Yvdot;
    m23 = -Yrdot;
    m33 = Iz - Nrdot;

    u = nu(1);
    v = nu(2);
    r = nu(3);

    % Row sums of C(nu)*nu + D(nu)*nu, matching tugboat3d.m's C and D
    f1 = -(m22*v + m23*r)*r + (-Xu - Xuu*abs(u))*u;
    f2 = m11*u*r + (-Yv - Yvv*abs(v))*v - Yr*r;
    f3 = (m22*v + m23*r)*u - m11*u*v - Nv*v - (Nr + Nrr*abs(r))*r;
end
