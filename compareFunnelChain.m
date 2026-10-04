%% COMPARE DURMAZ vs KINEMATIC CBF vs DYNAMIC HOCBF -- FUNNEL CHAIN
%
% Runs RSC.m to build a circular-funnel chain from start to goal (edit
% RSC.m's scenarioId / P.safetyMargin / P.minCircArea if you want a
% different map or a tighter chain), then drives the SAME chain with
% three controllers, all using Durmaz2024's own priority funnel-selection
% rule (highest-priority funnel -- i.e. closest to the goal -- containing
% the vessel):
%
%   Durmaz: plain Durmaz2024 law -> PI loops -> F (no safety filter; baseline)
%   CBF   : Durmaz2024 (u,w) -> kinematic CBF/HOCBF filter (CBF.m) ->
%           PI loops (lowLevelControl.m) -> F
%   HOCBF : Durmaz2024 (u,w) -> PI loops -> dynamic HOCBF force filter
%           (HOCBF.m) -> F
%
% Vessel: tugboat3d (default), otter3d or cybership3d. A caller can set a
% struct cfg before running this script to override vessel, U_m, profile
% and scaleR (see runShipComparison.m).

if ~exist('cfg', 'var'), cfg = struct(); end
clearvars -except cfg; clc; close all;

%% ---------------- Build the funnel chain (RSC.m, unmodified) ----------------
set(0, 'DefaultFigureVisible', 'off');
run('RSC.m');
set(0, 'DefaultFigureVisible', 'on');
close all;

if isempty(pathIds)
    error('RSC produced no path -- nothing to compare controllers on.');
end
fprintf('Funnel chain built: %d funnels, %d-funnel path.\n', numel(nodes), numel(pathIds));

%% ---------------- Simulation setup ----------------
dt       = 0.01;
goal_tol = 2.0;
vesselName = 'cybership';      % 'tugboat', 'otter' or 'cybership'
if isfield(cfg, 'vessel'), vesselName = cfg.vessel; end
model    = str2func([vesselName '3d']);   % xdot = model(x, [F_L; F_R])
vessel   = model();          % parameters and thrust limits of the vessel model
Flim     = [vessel.Fmin, vessel.Fmax];   % thrust limits per thruster [N]

%% ---------------- Mission speed (all controllers) ----------------
% Surge reference u = s(rho)*cos(alpha) with the speed magnitude s(rho):
%   MS.enable = false : s = k*rho                     (original law, k = 2*Kv)
%   MS.enable = true  : intermediate funnels  s = U_m (constant cruise speed)
%                       goal funnel           'linear': s = min(U_m, U_m*rho/R)
%                                             'tanh'  : s = U_m*tanh(k*rho/U_m)
% Any s(rho) >= 0 keeps the kinematic guarantees rho_dot = -s*cos(alpha)^2 <= 0
% and alpha_dot = -Ka*alpha.
% With MS.scaleR, U_m is replaced in every funnel by U_k = min(U_m, R_k/T_R):
% the heading converges in ~1/Ka s, during which the vessel travels U/Ka m,
% so T_R = 1/Ka keeps that turning distance within the funnel radius.
MS.enable  = true;
MS.U_m     = 5.0;            % mission (cruise) speed [m/s]
MS.profile = 'linear';       % goal-funnel slowdown: 'linear' or 'tanh'
MS.scaleR  = true;           % scale the cruise speed with the active funnel radius
for f = {'U_m', 'profile', 'scaleR'}
    if isfield(cfg, f{1}), MS.(f{1}) = cfg.(f{1}); end
end

%% ---------------- Durmaz2024 nominal law (CBF & HOCBF branches) ----------------
Kv = 0.05; Ka = 0.30; rho_tol = 0.05;
MS.T_R = 1/Ka;               % [s], U_k = min(U_m, R_k/T_R)

if MS.enable
    Rpath = [nodes(pathIds).radius];
    Upath = MS.U_m*ones(size(Rpath));
    if MS.scaleR, Upath = min(MS.U_m, Rpath/MS.T_R); end
    fprintf('Path funnels (start -> goal):\n  R [m]   = %s\n  U [m/s] = %s\n', ...
        mat2str(Rpath, 3), mat2str(Upath, 2));
end

% Simulation time: 4x the path length at the mission speed (at least 500 s)
Cpath = reshape([nodes(pathIds).c], 2, []);
Lpath = sum(vecnorm(diff(Cpath, 1, 2)));
Tsim  = 500;
if MS.enable, Tsim = max(Tsim, 4*Lpath/MS.U_m); end

%% ---------------- Path tracking (option) ----------------
% TRK.enable = false : Durmaz2024 steers to the active funnel's center.
% TRK.enable = true  : a path is generated from the RSC chain,
%   q_start -> c_1 -> c_2 -> ... -> c_goal (path funnel centers),
% and tracked with LOS guidance. Every segment c_k -> c_k+1 lies inside
% funnel k+1 (RSC puts each center inside its parent funnel), so the path
% stays inside the chain. The speed magnitude is the same s(rho) as above
% (mission speed, radius scaling, goal profile); only the steering changes.
TRK.enable    = false;
TRK.lookahead = 2.0;         % LOS lookahead distance [m]
TRK.Kpsi      = 0.5;         % heading gain: w_ref = Kpsi*(psi_LOS - psi) [1/s]
if isfield(cfg, 'track'), TRK.enable = cfg.track; end
TRK.wp = [q_start(:), Cpath];

%% ---------------- Disturbances ----------------
% None of them is known to the Durmaz law or to the CBF/HOCBF filters.
%   external : constant water current, added to the ground velocity
%   INS      : white noise on the measured position, heading and velocities;
%              funnel selection, controller and filters all use the measured
%              state, the margin h is logged from the true state
%   actuator : thrust efficiency of each thruster (e.g. a weak left motor)
%              and white noise on the applied thrust
DIST.Vc        = [0; 0];     % [m/s] current, world frame
DIST.ins_pos   = 0;          % [m] position noise std
DIST.ins_psi   = 0;          % [rad] heading noise std
DIST.ins_vel   = 0;          % [m/s, rad/s] noise std on u, v, r
DIST.act_gain  = [1; 1];     % applied/commanded thrust, [left; right]
DIST.act_noise = 0;          % [N] thrust noise std
DIST.seed      = 1;          % same noise sequence for every controller
if isfield(cfg, 'dist')
    for f = fieldnames(cfg.dist)'
        DIST.(f{1}) = cfg.dist.(f{1});
    end
end

%% ---------------- Low-level PI (lowLevelControl.m gains) ----------------
% Tugboat gains, scaled with the vessel's surge mass and yaw inertia so the
% loop bandwidths are the same on every vessel (factor 1 on the tugboat)
s_u = vessel.M(1,1)/11.07818;
s_r = vessel.M(3,3)/0.68525;
PI.P_speed = 100*s_u; PI.I_speed = 50*s_u;
PI.P_yaw   = 5*s_r;   PI.I_yaw   = 0.02*s_r;
PI.Ispeed_max = 2000*vessel.Fmax; PI.Iyaw_max = 20000;

%% ---------------- Kinematic CBF/HOCBF filter on (u,w) refs (CBF.m) ----------------
KF.k1 = 5; KF.k2 = 5;
KF.u_lim = max(1.0, MS.enable*MS.U_m); KF.w_lim = pi/2;   % must not undercut the mission speed
KF.use_omega_hocbf = true;
KF.slack_w = 1e4;
KF.qp_options = optimoptions('quadprog', 'Display', 'off');

%% ---------------- Dynamic HOCBF force filter (HOCBF.m) ----------------
HF.k1 = 5; HF.k2 = 5;
HF.qp_options = optimoptions('quadprog', 'Display', 'off');

%% ---------------- Run the three controllers on the same chain ----------------
ctrls = {'Durmaz', 'CBF', 'HOCBF'};
res = cell(1, numel(ctrls));
for c = 1:numel(ctrls)
    res{c} = runChainSim(ctrls{c}, nodes, pathIds, q_start, q_goal, dt, Tsim, ...
        goal_tol, Flim, Kv, Ka, rho_tol, PI, KF, HF, MS, DIST, model, vessel, TRK);
end

fprintf('\nVessel: %s', vessel.name);
if MS.enable
    fprintf('\nResults (mission speed U_m = %.2f m/s, goal profile = %s):\n', MS.U_m, MS.profile);
else
    fprintf('\nResults (original speed law s = k*rho):\n');
end
fprintf('%-6s | %-18s | %8s | %10s | %10s | %9s | %12s | %8s\n', ...
    'Ctrl', 'Outcome', 'T [s]', 'max viol', 't_out [s]', 'path [m]', 'int F^2 dt', 'filter %');
for c = 1:numel(ctrls)
    r = res{c};
    fprintf('%-6s | %-18s | %8.1f | %10.3f | %10.1f | %9.2f | %12.3g | %8.2f\n', ctrls{c}, ...
        tf2str(r.reached), r.time(end), r.viol, r.t_out, r.path_len, r.effort, 100*r.filt_frac);
end

%% ---------------- Plot 1: trajectories on the funnel chain ----------------
figure('Name','Funnel chain - trajectories','Color','w');
ax = gca; hold(ax, 'on'); axis(ax, 'equal');
xlim(ax, W(1:2)); ylim(ax, W(3:4));
for i = 1:numel(obs)
    plot(ax, obs{i}, 'FaceColor',[0 0 0], 'FaceAlpha',0.6, 'EdgeColor','none');
end
for k = 1:numel(pathIds)
    plot(ax, nodes(pathIds(k)).poly, 'FaceColor',[1 0.85 0.7], 'FaceAlpha',0.30, ...
        'EdgeColor',[1 0.5 0], 'LineWidth',1.0);
end
cols = lines(numel(ctrls));
lw = [4 2.5 1.2];                  % nested widths so identical paths stay visible
ls = {'-','--','-'};
h = gobjects(1, numel(ctrls));
for c = 1:numel(ctrls)
    h(c) = plot(ax, res{c}.X, res{c}.Y, 'Color', cols(c,:), 'LineWidth', lw(c), 'LineStyle', ls{c});
end
if TRK.enable
    plot(ax, TRK.wp(1,:), TRK.wp(2,:), 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
end
plot(ax, q_start(1), q_start(2), 'go', 'MarkerFaceColor','g');
plot(ax, q_goal(1),  q_goal(2),  'ro', 'MarkerFaceColor','r');
legend(h, ctrls, 'Location', 'bestoutside');
title(sprintf('%s: trajectory comparison', vessel.name));

%% ---------------- Plot 2: margin to the active funnel's boundary ----------------
figure('Name','Funnel chain - active funnel margin','Color','w');
hold on; grid on;
for c = 1:numel(ctrls)
    plot(res{c}.time, res{c}.h, 'Color', cols(c,:), 'LineWidth', lw(c), 'LineStyle', ls{c});
end
yline(0, 'k--');
xlabel('t [s]'); ylabel('h = R_{active} - \rho  [m]');
title(sprintf('%s: margin to the active funnel boundary (h < 0: outside)', vessel.name));
legend(ctrls, 'Location', 'best');

%% ---------------- Plot 3: surge speed u and yaw rate w ----------------
figure('Name','Funnel chain - surge speed and yaw rate','Color','w');
tiledlayout(2, 1, 'Padding','compact', 'TileSpacing','compact');
nexttile; hold on; grid on;
for c = 1:numel(ctrls)
    plot(res{c}.time, res{c}.u, 'Color', cols(c,:), 'LineWidth', lw(c), 'LineStyle', ls{c});
end
if MS.enable, yline(MS.U_m, 'k:'); end
ylabel('u [m/s]');
title(sprintf('%s: surge speed', vessel.name));
legend(ctrls, 'Location', 'best');
nexttile; hold on; grid on;
for c = 1:numel(ctrls)
    plot(res{c}.time, res{c}.w, 'Color', cols(c,:), 'LineWidth', lw(c), 'LineStyle', ls{c});
end
xlabel('t [s]'); ylabel('w [rad/s]');
title('Yaw rate');

%% ================================================================
%% Local functions
%% ================================================================

function o = runChainSim(ctrl, nodes, pathIds, q_start, q_goal, dt, Tsim, ...
        goal_tol, Flim, Kv, Ka, rho_tol, PI, KF, HF, MS, DIST, model, vessel, TRK)

    N = round(Tsim/dt) + 1;
    time = (0:N-1)*dt;

    % Initial state array [X Y psi u v r F_L F_R], starts at rest at q_start
    x = [q_start(1); q_start(2); 0; 0; 0; 0; 0; 0];

    % Path funnels (for the margin when the vessel is outside all of them)
    Cp = reshape([nodes(pathIds).c], 2, []);
    Rp = [nodes(pathIds).radius];

    % Logs
    X_h = zeros(1,N); Y_h = zeros(1,N); h_h = zeros(1,N);
    u_h = zeros(1,N); w_h = zeros(1,N);
    F_h = zeros(2,N); filt_h = false(1,N);
    rng(DIST.seed);
    int_speed = 0; int_yaw = 0;
    last = N;
    reached = false;
    seg = 1;                                    % current path segment (TRK)

    % Time loop
    for k = 1:N
        % Measured state (INS noise); the controller only sees xm
        xm = x;
        xm(1:2) = x(1:2) + DIST.ins_pos*randn(2,1);
        xm(3)   = wrapToPiLocal(x(3) + DIST.ins_psi*randn);
        xm(4:6) = x(4:6) + DIST.ins_vel*randn(3,1);

        % Active funnel: highest-priority funnel containing the vessel
        position = xm(1:2)';
        [~, active_node_id] = selectActiveFunnel(position, nodes, pathIds);
        ctr = nodes(active_node_id).c(:);
        R   = nodes(active_node_id).radius;
        is_goal = (active_node_id == pathIds(end));
        [rho_m, ~] = polarState(xm, ctr);

        % True distance to the active funnel center and margin h = R - rho (log it)
        [rho, ~] = polarState(x, ctr);
        X_h(k) = x(1); Y_h(k) = x(2); h_h(k) = R - rho;
        out_all = vecnorm(Cp - x(1:2)) - Rp;        % > 0: outside path funnel j
        if all(out_all > 0), h_h(k) = -min(out_all); end   % outside all: distance to the nearest
        u_h(k) = x(4); w_h(k) = x(6);                      % surge speed, yaw rate

        % Stop when the goal is reached
        if norm(x(1:2)' - q_goal) <= goal_tol
            reached = true; last = k; break;
        end
        if k == N, break; end

        % Thrust F = [F_L; F_R]: Durmaz2024 nominal law + PI (+ filter)
        if TRK.enable
            s_mag = speedProfile(rho_m, R, 2*Kv, MS, is_goal);
            [un, wn, seg] = losNominal(xm, TRK, seg, s_mag);
        else
            [un, wn] = durmazNominal(xm, ctr, R, Kv, Ka, rho_tol, MS, is_goal);
        end
        u_ref = un; w_ref = wn;
        % CBF: correct the references before they reach the PI loops
        if strcmp(ctrl, 'CBF') && rho_m > rho_tol
            [u_ref, w_ref] = kinFilter(xm, un, wn, ctr, R, KF);
            filt_h(k) = abs(u_ref - un) > 1e-4 || abs(w_ref - wn) > 1e-4;
        end

        % PI loops on surge and yaw-rate errors
        e_speed = u_ref - xm(4);
        e_yaw   = w_ref - xm(6);
        int_speed_new = min(max(int_speed + e_speed*dt, -PI.Ispeed_max), PI.Ispeed_max);
        int_yaw_new   = min(max(int_yaw   + e_yaw  *dt, -PI.Iyaw_max),   PI.Iyaw_max);

        tauX = PI.P_speed*e_speed + PI.I_speed*int_speed_new;
        tauN = PI.P_yaw*e_yaw + PI.I_yaw*int_yaw_new;

        % Thrust allocation: (tauX, tauN) -> (F_L, F_R), d = moment arm
        d = vessel.d;
        FLc = tauX/2 + tauN/(2*d);
        FRc = tauX/2 - tauN/(2*d);
        FL = min(max(FLc, Flim(1)), Flim(2));
        FR = min(max(FRc, Flim(1)), Flim(2));

        % Anti-windup (conditional integration): keep the integrators frozen
        % while a thruster saturates. (Resetting them so that the output
        % equals the saturated thrust drove int_yaw to a large opposite bias
        % whenever the P term alone exceeded the thrust limits; with
        % I_yaw = 0.02 it took minutes to unwind and blocked the yaw loop.)
        if (FL == FLc) && (FR == FRc)
            int_speed = int_speed_new;
            int_yaw   = int_yaw_new;
        end

        F = [FL; FR];
        % HOCBF: correct the thrust itself
        if strcmp(ctrl, 'HOCBF') && rho_m > rho_tol
            F_nom = F;
            F = hocbfForceFilter(xm, F, ctr, R, HF, Flim, model, vessel);
            filt_h(k) = norm(F - F_nom) > 1e-3;
        end

        F_h(:,k) = F;

        % Actuator disturbance: thrust efficiency and noise
        Fa = DIST.act_gain.*F + DIST.act_noise*randn(2,1);

        % Plant step (RK4, thrust held constant over dt)
        cur = [DIST.Vc; zeros(6,1)];                    % current moves X, Y only
        f1 = model(x, Fa) + cur;
        x2 = x + 0.5*dt*f1;
        f2 = model(x2, Fa) + cur;
        x3 = x + 0.5*dt*f2;
        f3 = model(x3, Fa) + cur;
        x4 = x + dt*f3;
        f4 = model(x4, Fa) + cur;
        x  = x + (dt/6)*(f1 + 2*f2 + 2*f3 + f4);
        x(3) = wrapToPiLocal(x(3));
    end

    % Results: trim logs; viol = worst excursion outside the active funnel
    idx = 1:last;
    o.time = time(idx); o.X = X_h(idx); o.Y = Y_h(idx); o.h = h_h(idx);
    o.u = u_h(idx); o.w = w_h(idx);
    o.filt_frac = mean(filt_h(idx));                       % share of steps the filter changed the command
    o.reached = reached;
    o.viol = max(0, -min(o.h));
    o.t_out    = dt*sum(o.h < 0);                          % time spent outside the active funnel
    o.path_len = sum(hypot(diff(o.X), diff(o.Y)));
    o.effort   = dt*sum(sum(F_h(:,idx).^2));               % integral of F_L^2 + F_R^2
end

function [active_path_idx, active_node_id] = selectActiveFunnel(position, nodes, pathIds)
% Durmaz2024.m's rule: the highest-priority (goal-most) funnel whose disc
% contains the vessel; pathIds is ordered start-side -> goal-side.
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

function [rho, alpha] = polarState(x, ctr)
    ex = x(1) - ctr(1);
    ey = x(2) - ctr(2);
    rho = hypot(ex, ey);
    phi = atan2(-ey, -ex);
    alpha = wrapToPiLocal(phi - x(3));
end

function [un, wn] = durmazNominal(x, ctr, R, Kv, Ka, rho_tol, MS, is_goal)
% Durmaz2024 Eq. 33 with the surge gain 2*Kv*rho replaced by s(rho)
% (speedProfile); wn keeps its form so alpha_dot = -Ka*alpha still holds.
    [rho, alpha] = polarState(x, ctr);
    if rho <= rho_tol
        un = 0; wn = 0;
    else
        un = speedProfile(rho, R, 2*Kv, MS, is_goal)*cos(alpha);
        wn = Ka*alpha + (un/rho)*sin(alpha);
    end
end

function [un, wn, seg] = losNominal(x, TRK, seg, s)
% LOS path following on the polyline TRK.wp (2 x n). The segment index only
% moves forward: it advances once the vessel's projection passes the end of
% the segment. Desired heading psi_LOS = chi_p - atan(e/lookahead), with e the
% cross-track error (positive to the left of the path).
    p = x(1:2);
    nseg = size(TRK.wp, 2) - 1;
    while seg < nseg
        a = TRK.wp(:, seg); t = TRK.wp(:, seg+1) - a;
        if dot(t, t) < 1e-12 || dot(p - a, t)/dot(t, t) >= 1
            seg = seg + 1;
        else
            break;
        end
    end
    a = TRK.wp(:, seg); b = TRK.wp(:, seg+1);
    Lseg = norm(b - a); t = (b - a)/max(Lseg, 1e-12);
    if seg == nseg && Lseg - dot(p - a, t) < TRK.lookahead
        psi_d = atan2(b(2) - p(2), b(1) - p(1));    % end of the path: aim at the goal
    else
        chi_p = atan2(t(2), t(1));
        e = t(1)*(p(2) - a(2)) - t(2)*(p(1) - a(1));
        psi_d = chi_p - atan(e/TRK.lookahead);
    end
    err = wrapToPiLocal(psi_d - x(3));
    un = s*max(cos(err), 0);                    % no reversing, slow down while turning
    wn = TRK.Kpsi*err;
end

function s = speedProfile(rho, R, k, MS, is_goal)
% Surge speed magnitude s(rho) (see MS above).
    U = MS.U_m;
    if MS.scaleR, U = min(U, R/MS.T_R); end         % cruise speed of this funnel
    if ~MS.enable                                   % original law
        s = k*rho;
    elseif ~is_goal                                 % cruise through the chain
        s = U;
    elseif strcmp(MS.profile, 'linear')             % Durmaz law with 2*Kv = U/R
        s = min(U, U*rho/R);
    else                                            % 'tanh'
        s = U*tanh(k*rho/U);
    end
end

function [u_ref, w_ref] = kinFilter(x, u_nom, w_nom, ctr, R, KF)
% Kinematic CBF / HOCBF filter on the references (u_ref, w_ref).
%
% Model (state xi = [X; Y; psi], input nu = [u; w]):
%     xi_dot = f(xi) + g(xi)*nu
%     f = R(psi)*[0; v]        (measured sway = drift, not controllable)
%     g = [cos psi 0; sin psi 0; 0 1]
% Barrier:  b(xi) = R - rho   (b >= 0 inside the funnel)
%
%   Lf b = v sin(alpha)            Lg b = [cos(alpha), 0]        -> rel. degree 1 in u, 2 in w
%
% (1) CBF on u:      Lf b + Lg b*nu + k1*b >= 0
% (2) HOCBF on w:    Lf2 b + LgLf b*nu + (k1+k2)*Lf b + k1*k2*b >= 0
%     with surge frozen at its measured value (drift ftilde = f + g_u*u_meas):
%         Lftilde b = u cos(alpha) + v sin(alpha)   (= bdot)
%         Lftilde^2 b = -w_t^2/rho      LgLftilde b = [0, -w_t]
%     (u_dot, v_dot neglected: kinematic-level model)
% QP:  min ||nu - nu_nom||^2 + Wd*(dA^2 + dB^2)  s.t. (1),(2) with slacks, box bounds.

    [rho, alpha] = polarState(x, ctr);
    u = x(4); v = x(5);
    wt = -u*sin(alpha) + v*cos(alpha);          % tangential velocity

    b = R - rho;

    % --- first-order Lie derivatives
    Lfb = v*sin(alpha);
    Lgb = [cos(alpha), 0];

    % --- second-order Lie derivatives (surge quasi-static)
    Lftb  = Lfb + Lgb(1)*u;                     % = bdot
    Lf2b  = -wt^2/rho;
    LgLfb = [0, -wt];

    % --- constraints in the form  A_i*nu <= c_i
    A1 = -Lgb;      c1 = Lfb + KF.k1*b;
    A2 = -LgLfb;    c2 = Lf2b + (KF.k1 + KF.k2)*Lftb + KF.k1*KF.k2*b;

    nu_nom = [u_nom; w_nom];
    viol1 = A1*nu_nom > c1;
    viol2 = KF.use_omega_hocbf && (A2*nu_nom > c2);

    u_ref = u_nom; w_ref = w_nom;
    if ~viol1 && ~viol2
        return;                                  % nominal input already safe
    end

    % --- CBF-QP with slacks, decision vector z = [u; w; d1; d2]
    Wd = KF.slack_w;
    H = diag([2, 2, 2*Wd, 2*Wd]);
    f = [-2*nu_nom; 0; 0];

    A = [A1, -1,  0];
    c = c1;
    if KF.use_omega_hocbf
        A = [A; A2, 0, -1];
        c = [c; c2];
    end

    lb = [-KF.u_lim; -KF.w_lim; 0;   0];
    ub = [ KF.u_lim;  KF.w_lim; inf; inf];

    [z, ~, exitflag] = quadprog(H, f, A, c, [], [], lb, ub, [nu_nom; 0; 0], KF.qp_options);
    if exitflag > 0
        u_ref = z(1);
        w_ref = z(2);
    end
end

function [b, Lfb, Lf2b, LgLfb] = funnelLie(x, ctr, R, model, vessel)
% Dynamic model (state xi = [X;Y;psi;u;v;r], input F = [F_L; F_R]):
%     xi_dot = f(xi) + g(xi)*F      (thrust enters through nu_dot = M^-1 (B F - (C+D) nu))
% Barrier:  b(xi) = R - rho
%
%   Lg b = 0                          (force does not appear in b_dot)
%   Lf b = u cos(a) + v sin(a)        (= b_dot)
%   ->  relative degree 2 in F:
%   b_ddot = Lf2 b + LgLf b * F
%     Lf2 b   = cos(a)*d1 + sin(a)*d2 - w_t^2/rho - w_t*r ,  d = -M^-1 (C nu + D nu)
%     LgLf b  = cos(a)*G(1,:) + sin(a)*G(2,:),               G = M^-1 B
%   with a = alpha (heading error to the center), w_t = -u sin(a) + v cos(a).
% M, C(nu), D(nu) come from the vessel model: with zero thrust its nu_dot is
% exactly d = -M^-1 (C nu + D nu), and G = M^-1 B_prop from its parameters.
    nu = x(4:6);
    xd = model([0; 0; 0; nu; 0; 0], [0; 0]);
    d  = xd(4:6);                               % drift acceleration of nu
    G  = vessel.M \ vessel.B_prop;

    [rho, alpha] = polarState(x, ctr);
    u = nu(1); v = nu(2); r = nu(3);
    wt = -u*sin(alpha) + v*cos(alpha);

    b     = R - rho;
    Lfb   = u*cos(alpha) + v*sin(alpha);
    Lf2b  = cos(alpha)*d(1) + sin(alpha)*d(2) - wt^2/rho - wt*r;
    LgLfb = cos(alpha)*G(1,:) + sin(alpha)*G(2,:);          % 1x2 row, acts on F
end

function F = hocbfForceFilter(x, F_nom, ctr, R, HF, Flim, model, vessel)
% 2nd-order HOCBF-QP on thrust:
%     min ||F - F_nom||^2
%     s.t. Lf2 b + LgLf b*F + (k1+k2)*Lf b + k1*k2*b >= 0 ,   Flim(1) <= F <= Flim(2)
% i.e.  -LgLf b * F <= Lf2 b + (k1+k2) Lf b + k1 k2 b.
    [b, Lfb, Lf2b, LgLfb] = funnelLie(x, ctr, R, model, vessel);
    A = -LgLfb;
    c = Lf2b + (HF.k1 + HF.k2)*Lfb + HF.k1*HF.k2*b;

    F = F_nom;
    if A*F_nom > c                              % nominal thrust violates the constraint
        [Fq, ~, exitflag] = quadprog(2*eye(2), -2*F_nom, A, c, [], [], ...
            [Flim(1); Flim(1)], [Flim(2); Flim(2)], F_nom, HF.qp_options);
        if exitflag > 0
            F = Fq;
        end
    end
end

function a = wrapToPiLocal(a)
    a = mod(a + pi, 2*pi) - pi;
end

function s = tf2str(tf)
    if tf, s = 'reached goal'; else, s = 'did not reach goal'; end
end
