%% FEEDBACK LINEARIZATION + (ROBUST) CONTROL BARRIER FUNCTION SAFETY FILTER
%
% Durmaz2024.m's nominal circular-funnel navigation law (Eq. 33), UNCHANGED,
% wrapped in a per-step safety filter following Wang, Xiao, Gonzalez-Garcia,
% Swevers, Ratti & Rus, "Robust Model Predictive Control with Control
% Barrier Functions for Autonomous Surface Vessels," ICRA 2024 -- their
% Def. 1-3 / Theorem 1 (HOCBF, from Xiao & Belta) and robust-CBF
% extension (their Eq. 20-24), applied here to circular funnel-boundary
% safety instead of their tube-around-a-trajectory application, and fused
% with feedback linearization instead of their MPC.
%
% Relative degree: for the KINEMATIC unicycle (this file) the barrier
%   b(x) = R^2 - rho^2,  rho = ||position - funnel_center||
% has relative degree 1 in v (differentiate once, v appears directly --
% see derivation below), so this is the m=1 case of their Def. 3, which
% reduces to an ordinary CBF (their own remark after Def. 3). The DYNAMIC
% tugboat version needs true relative-degree-2 HOCBF, since force is two
% integrations from position; that is a separate file.
%
% Derivation (control-affine unicycle, x=[X;Y;psi], u=[v;omega], plus a
% bounded external disturbance w = [wx;wy;0] with |w| <= wbar):
%
%   xdot = f(x) + g(x)*u + w,   f(x)=0,   g(x) = [cos(psi) 0; sin(psi) 0; 0 1]
%
%   b(x)    = R^2 - (X-cx)^2 - (Y-cy)^2
%   db/dx   = [-2(X-cx), -2(Y-cy), 0]
%   L_f b   = 0
%   L_g b   = [-2(X-cx)*cos(psi) - 2(Y-cy)*sin(psi), 0]   -- only couples to v
%
% Requiring L_f b + L_g b*u + alpha_1(b) >= 0 to hold for every admissible
% disturbance (worst case, Cauchy-Schwarz over |w|<=wbar) gives the ROBUST
% CBF constraint actually enforced here (their Eq. 21, alpha_1 taken
% linear, alpha_1(b) = gamma*b):
%
%   Lgb_v * v + gamma*b - 2*rho*wbar >= 0,      Lgb_v = L_g b (v-component)
%
% Since only v appears, the minimally-invasive safety filter
%   min ||u - u_nom||^2   s.t. the constraint above
% leaves omega untouched and only ever corrects v -- solved each step as
% a 2-variable QP via quadprog (matching the paper's own Eq. 24 QP
% formulation) rather than hand-deriving the 1D projection case-split,
% since a general QP solver is less error-prone than a hand-derived
% closed form (mirroring the deliberate use of exact symbolic/verified
% derivations elsewhere in this project).
%
% GAP THIS CLOSES (see compare_fbl_hocbf.m for the demonstration):
% Durmaz2024.m's rho_dot <= 0 proof assumes perfect, undisturbed
% kinematics. Under a constant disturbance of magnitude |w|, the nominal
% law's commanded restoring speed (~2*Kv*rho) balances the disturbance at
% a STEADY-STATE distance rho_ss ~ |w|/(2*Kv) -- a permanent funnel-
% boundary violation whenever |w| > 2*Kv*R, not a transient blip. The
% robust CBF term above is designed to counter exactly this: since the
% kinematic model has no actuator limit, the QP can always find enough v
% near the boundary to hold rho <= R despite the same disturbance.

%% GAINS (Durmaz2024.m's own nominal law, unchanged, plus the CBF's own)

if ~exist('Kv','var'), Kv = 0.10; end
if ~exist('Ka','var'), Ka = 1.0;  end

if ~exist('gamma_cbf','var'), gamma_cbf = 1.0; end % class-K gain alpha_1(b) = gamma_cbf*b
if ~exist('wbar','var'),      wbar      = 0.0; end % disturbance bound the CBF defends against [m/s]
if ~exist('cbf_enabled','var'), cbf_enabled = true; end % false => bypass the QP, pure nominal law

% Actual disturbance applied to the simulated plant (separate from wbar,
% the BOUND the CBF assumes -- set wbar >= norm([wx,wy]) for the robust
% guarantee to actually hold).
if ~exist('wx','var'), wx = 0.0; end
if ~exist('wy','var'), wy = 0.0; end
if ~exist('disturbance_radius','var'), disturbance_radius = inf; end
if ~exist('ramp_duration','var'), ramp_duration = 0; end
% ramp_duration: seconds over which the ACTUAL disturbance grows linearly
% from 0 to its full (wx,wy) magnitude after latching on, instead of
% jumping to full strength instantly (ramp_duration=0, the default).
% wbar (the CBF's assumed worst-case bound) is NOT ramped -- per the
% robust-CBF theory it represents a fixed, known bound the controller
% defends against from the start, which remains a valid (if temporarily
% conservative) bound throughout the ramp-up, since the actual
% disturbance never exceeds its own eventual full magnitude.
% (wx,wy) LATCH on permanently the first time the vehicle comes within
% disturbance_radius of the FINAL goal point, and never turn back off --
% models a current that's present in the final approach/passage, without
% the self-canceling feedback a REVERSIBLE spatial gate produces: pushing
% the vehicle back out of the zone would normally remove the disturbance,
% letting it drift back in, get pushed out again, and only ever hover at
% the gate boundary instead of showing a clean violation (tried first,
% same failure mode as gating by the ACTIVE funnel). A one-way latch has
% no such feedback. Default inf means "latch immediately" (t=0), i.e.
% uniform disturbance for the whole run -- which is itself fine for a
% single funnel, but for a multi-funnel CHAIN a uniform disturbance
% typically prevents the vehicle from ever leaving the FIRST funnel in
% the first place (any disturbance strong enough to violate a small
% final funnel also traps the vehicle in an orbit around whichever
% funnel is currently active, since progression to the next funnel
% requires reducing distance to ITS center against the same
% disturbance). Setting disturbance_radius finite lets the vehicle
% traverse the whole chain undisturbed and only latches the disturbance
% on once genuinely close to the end.

if ~exist('dt_sim','var'),   dt_sim   = 0.01;  end
if ~exist('sim_time','var'), sim_time = 200.0; end
if ~exist('goal_tol','var'), goal_tol = 0.05;  end
if ~exist('rho_arrival_tol','var'), rho_arrival_tol = 0.05; end

theta0 = 0.0;

qp_options = optimoptions('quadprog', 'Display', 'off');

%% INITIAL STATE

% state = [X; Y; psi]
state = [q_start(1); q_start(2); theta0];

time = 0:dt_sim:sim_time;
N = numel(time);

state_hist          = zeros(N, 3);
v_hist               = zeros(N, 1);
omega_hist           = zeros(N, 1);
vnom_hist            = zeros(N, 1); % nominal (pre-QP) surge command, for comparison
rho_hist             = zeros(N, 1);
b_hist               = zeros(N, 1); % barrier value R^2 - rho^2 (>=0 means safe)
active_funnel_hist   = zeros(N, 1);
funnel_violation     = false(N, 1);
overshoot_hist       = zeros(N, 1);
cbf_active_hist      = false(N, 1); % true whenever the QP actually altered v

state_hist(1,:) = state.';

%% SIMULATION

last_idx = N;
disturbance_latched = false; % one-way: set true once, never reset
latch_time = NaN; % simulation time at which the latch fired, for ramping

for k = 1:N-1

    position = state(1:2).';
    psi = state(3);

    % Select highest-priority funnel containing current vehicle position.
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

    % One-way latch: once within disturbance_radius of the goal POINT, OR
    % once the goal FUNNEL itself first becomes active (whichever happens
    % first), the disturbance stays on for the rest of the run (see the
    % header note above). The funnel-based trigger is the more robust of
    % the two on a multi-funnel chain: a raw distance threshold can fire
    % while an earlier, non-final funnel is still active (if that funnel
    % happens to extend within the threshold distance of the goal point),
    % trapping the vehicle there instead of the intended final funnel.
    if norm(state(1:2).' - q_goal) <= disturbance_radius || is_goal_funnel
        if ~disturbance_latched
            latch_time = time(k);
        end
        disturbance_latched = true;
    end

    % Bearing toward funnel center, Durmaz2024.m's Eq. (11)/(17) convention
    phi = atan2(-y, -x);
    alpha = wrapToPiLocal(phi - psi);

    % Nominal feedback-linearizing law, Durmaz2024.m's Eq. (33), UNCHANGED
    if is_goal_funnel && rho <= rho_arrival_tol
        v_nom = 0;
        omega_nom = 0;
    else
        v_nom     = 2 * Kv * rho * cos(alpha);
        omega_nom = Ka * alpha + (v_nom / rho) * sin(alpha);
    end

    % Robust CBF safety filter: min ||u - u_nom||^2 s.t. the constraint.
    % cbf_enabled=false bypasses the QP entirely (pure nominal law) --
    % used to ablate the filter's effect with everything else identical.
    if cbf_enabled
        Lgb_v = -2*x*cos(psi) - 2*y*sin(psi);

        if disturbance_latched
            wbar_active = wbar;
        else
            wbar_active = 0; % nothing to be robust against before the disturbance latches on
        end

        H = 2*eye(2);
        f = -2*[v_nom; omega_nom];
        A = [-Lgb_v, 0];
        bineq = gamma_cbf*b - 2*rho*wbar_active;

        [u_opt, ~, exitflag] = quadprog(H, f, A, bineq, [], [], [], [], [v_nom; omega_nom], qp_options);

        if exitflag > 0
            v = u_opt(1);
            omega = u_opt(2);
        else
            % QP infeasible (can happen if already unsafe and Lgb_v == 0):
            % fall back to the nominal command rather than fail silently.
            v = v_nom;
            omega = omega_nom;
        end
    else
        v = v_nom;
        omega = omega_nom;
    end

    cbf_active_hist(k) = abs(v - v_nom) > 1e-9;

    % Disturbed kinematics: nominal unicycle + external drift, ramping
    % linearly from 0 to full (wx,wy) strength over ramp_duration seconds
    % after latching on (ramp_duration=0 -> instant full strength).
    if disturbance_latched
        if ramp_duration > 0
            ramp_scale = min((time(k) - latch_time) / ramp_duration, 1);
        else
            ramp_scale = 1;
        end
        wx_active = ramp_scale * wx;
        wy_active = ramp_scale * wy;
    else
        wx_active = 0; wy_active = 0;
    end
    xdot = v*cos(psi) + wx_active;
    ydot = v*sin(psi) + wy_active;

    state = state + dt_sim*[xdot; ydot; omega];
    state(3) = wrapToPiLocal(state(3));

    % Save
    state_hist(k+1,:)         = state.';
    v_hist(k)                 = v;
    omega_hist(k)               = omega;
    vnom_hist(k)                  = v_nom;
    rho_hist(k)                     = rho;
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
vnom_hist = vnom_hist(1:last_idx);
rho_hist = rho_hist(1:last_idx);
b_hist = b_hist(1:last_idx);
active_funnel_hist = active_funnel_hist(1:last_idx);
funnel_violation = funnel_violation(1:last_idx);
overshoot_hist = overshoot_hist(1:last_idx);
cbf_active_hist = cbf_active_hist(1:last_idx);

if last_idx > 1
    v_hist(end) = v_hist(end-1);
    omega_hist(end) = omega_hist(end-1);
    vnom_hist(end) = vnom_hist(end-1);
    rho_hist(end) = rho_hist(end-1);
    b_hist(end) = b_hist(end-1);
    active_funnel_hist(end) = active_funnel_hist(end-1);
end

fprintf('\n--- Feedback Linearization + Robust CBF Safety Filter ---\n');
fprintf('Kv = %.3f, Ka = %.3f, gamma_cbf = %.3f, wbar = %.3f, disturbance=[%.3f,%.3f]\n', ...
    Kv, Ka, gamma_cbf, wbar, wx, wy);
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

title(ax1, 'Funnel Chain + Closed-Loop Trajectory (FBL + Robust CBF)');

%% PLOT 2: BARRIER VALUE AND CBF ACTIVATION

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

%% PLOT 3: CONTROL INPUTS

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
