%% COMPARE PFL vs KINEMATIC CBF vs DYNAMIC HOCBF -- FUNNEL CHAIN
%
% Runs RSC.m to build a circular-funnel chain from start to goal (edit
% RSC.m's scenarioId / P.safetyMargin / P.minCircArea if you want a
% different map or a tighter chain), then drives the SAME chain with
% three controllers, all using Durmaz2024's own priority funnel-selection
% rule (highest-priority funnel -- i.e. closest to the goal -- containing
% the vessel):
%
%   Durmaz: plain Durmaz2024 law -> PI loops -> F (no safety filter; baseline)
%   PFL   : SwayAwarePFL.m's own exact partial-feedback-linearization
%           funnel-tracking law (no filter, own funnel controller)
%   CBF   : Durmaz2024 (u,w) -> kinematic CBF/HOCBF filter (CBF.m) ->
%           PI loops (lowLevelControl.m) -> F
%   HOCBF : Durmaz2024 (u,w) -> PI loops -> dynamic HOCBF force filter
%           (HOCBF.m) -> F

clear; clc; close all;

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
Tsim     = 500;
goal_tol = 2.0;
Fmax     = 1000;             % thruster force limit [N], shared by all three

%% ---------------- Durmaz2024 nominal law (CBF & HOCBF branches) ----------------
Kv = 0.05; Ka = 0.30; rho_tol = 0.05;

%% ---------------- Low-level PI (lowLevelControl.m gains) ----------------
PI.P_speed = 100; PI.I_speed = 50;
PI.P_yaw   = 5;   PI.I_yaw   = 0.02;
PI.Ispeed_max = 2000*Fmax; PI.Iyaw_max = 20000;

%% ---------------- Kinematic CBF/HOCBF filter on (u,w) refs (CBF.m) ----------------
KF.k1 = 5; KF.k2 = 5;
KF.u_lim = 1.0; KF.w_lim = pi/2;
KF.use_omega_hocbf = true;
KF.slack_w = 1e4;
KF.qp_options = optimoptions('quadprog', 'Display', 'off');

%% ---------------- Dynamic HOCBF force filter (HOCBF.m) ----------------
HF.k1 = 5; HF.k2 = 5;
HF.qp_options = optimoptions('quadprog', 'Display', 'off');

%% ---------------- PFL gains (SwayAwarePFL.m defaults) ----------------
PFLg.k_rho = 0.10; PFLg.k_alpha = 0.30;
PFLg.k_u   = 1.0;  PFLg.k_r     = 1.0;
PFLg.rho_floor = 1.0;

%% ---------------- Run the three controllers on the same chain ----------------
ctrls = {'Durmaz', 'PFL', 'CBF', 'HOCBF'};
res = cell(1, numel(ctrls));
for c = 1:numel(ctrls)
    res{c} = runChainSim(ctrls{c}, nodes, pathIds, q_start, q_goal, dt, Tsim, ...
        goal_tol, Fmax, Kv, Ka, rho_tol, PI, KF, HF, PFLg);
end

fprintf('\nResults:\n');
fprintf('%-6s | %-18s | %8s | %10s | %10s | %9s | %12s\n', ...
    'Ctrl', 'Outcome', 'T [s]', 'max viol', 't_out [s]', 'path [m]', 'int F^2 dt');
for c = 1:numel(ctrls)
    r = res{c};
    fprintf('%-6s | %-18s | %8.1f | %10.3f | %10.1f | %9.2f | %12.3g\n', ctrls{c}, ...
        tf2str(r.reached), r.time(end), r.viol, r.t_out, r.path_len, r.effort);
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
lw = [5 3.5 2.5 1.2];            % nested widths so identical paths stay visible
ls = {'-','-','--','-'};
h = gobjects(1, numel(ctrls));
for c = 1:numel(ctrls)
    h(c) = plot(ax, res{c}.X, res{c}.Y, 'Color', cols(c,:), 'LineWidth', lw(c), 'LineStyle', ls{c});
end
plot(ax, q_start(1), q_start(2), 'go', 'MarkerFaceColor','g');
plot(ax, q_goal(1),  q_goal(2),  'ro', 'MarkerFaceColor','r');
legend(h, ctrls, 'Location', 'bestoutside');
title('Funnel chain: trajectory comparison');

%% ---------------- Plot 2: margin to the active funnel's boundary ----------------
figure('Name','Funnel chain - active funnel margin','Color','w');
hold on; grid on;
for c = 1:numel(ctrls)
    plot(res{c}.time, res{c}.h, 'Color', cols(c,:), 'LineWidth', lw(c), 'LineStyle', ls{c});
end
yline(0, 'k--');
xlabel('t [s]'); ylabel('h = R_{active} - \rho  [m]');
title('Margin to the active funnel boundary (h < 0: outside)');
legend(ctrls, 'Location', 'best');

%% ================================================================
%% Local functions
%% ================================================================

function o = runChainSim(ctrl, nodes, pathIds, q_start, q_goal, dt, Tsim, ...
        goal_tol, Fmax, Kv, Ka, rho_tol, PI, KF, HF, PFLg)

    N = round(Tsim/dt) + 1;
    time = (0:N-1)*dt;

    % Initial state array [X Y psi u v r F_L F_R], starts at rest at q_start
    x = [q_start(1); q_start(2); 0; 0; 0; 0; 0; 0];

    % Logs
    X_h = zeros(1,N); Y_h = zeros(1,N); h_h = zeros(1,N);
    F_h = zeros(2,N);
    int_speed = 0; int_yaw = 0;
    last = N;
    reached = false;

    % Time loop
    for k = 1:N
        % Active funnel: highest-priority funnel containing the vessel
        position = x(1:2)';
        [~, active_node_id] = selectActiveFunnel(position, nodes, pathIds);
        ctr = nodes(active_node_id).c(:);
        R   = nodes(active_node_id).radius;

        % Distance to active funnel center and margin h = R - rho (log it)
        [rho, ~] = polarState(x, ctr);
        X_h(k) = x(1); Y_h(k) = x(2); h_h(k) = R - rho;

        % Stop when the goal is reached
        if norm(x(1:2)' - q_goal) <= goal_tol
            reached = true; last = k; break;
        end
        if k == N, break; end

        % Compute thrust F = [F_L; F_R] according to the controller
        switch ctrl
            case 'PFL'
                F = pflForce(x, ctr, PFLg, Fmax);

            otherwise % 'CBF' or 'HOCBF': Durmaz2024 nominal law + PI
                % Durmaz2024 nominal surge / yaw-rate references
                [un, wn] = durmazNominal(x, ctr, Kv, Ka, rho_tol);
                u_ref = un; w_ref = wn;
                % CBF: correct the references before they reach the PI loops
                if strcmp(ctrl, 'CBF') && rho > rho_tol
                    [u_ref, w_ref] = kinFilter(x, un, wn, ctr, R, KF);
                end

                % PI loops on surge and yaw-rate errors
                e_speed = u_ref - x(4);
                e_yaw   = w_ref - x(6);
                int_speed = min(max(int_speed + e_speed*dt, -PI.Ispeed_max), PI.Ispeed_max);
                int_yaw   = min(max(int_yaw   + e_yaw  *dt, -PI.Iyaw_max),   PI.Iyaw_max);

                tauX = PI.P_speed*e_speed + PI.I_speed*int_speed;
                tauN = PI.P_yaw*e_yaw + PI.I_yaw*int_yaw;

                % Thrust allocation: (tauX, tauN) -> (F_L, F_R), d = moment arm
                d = 0.29/2;
                FLc = tauX/2 + tauN/(2*d);
                FRc = tauX/2 - tauN/(2*d);
                FL = min(max(FLc, -Fmax), Fmax);
                FR = min(max(FRc, -Fmax), Fmax);

                % Anti-windup: unwind the integrators if thrust saturated
                if (FL ~= FLc) || (FR ~= FRc)
                    tauX_ach = FL + FR;
                    tauN_ach = d*(FL - FR);
                    int_speed = int_speed + (tauX_ach - tauX)/PI.I_speed;
                    int_speed = min(max(int_speed, -PI.Ispeed_max), PI.Ispeed_max);
                    int_yaw = int_yaw + (tauN_ach - tauN)/PI.I_yaw;
                    int_yaw = min(max(int_yaw, -PI.Iyaw_max), PI.Iyaw_max);
                end

                F = [FL; FR];
                % HOCBF: correct the thrust itself
                if strcmp(ctrl, 'HOCBF') && rho > rho_tol
                    F = hocbfForceFilter(x, F, ctr, R, HF, Fmax);
                end
        end

        F_h(:,k) = F;

        % Plant step (RK4, thrust held constant over dt)
        f1 = tugboat3d(x, F);
        x2 = x + 0.5*dt*f1;
        f2 = tugboat3d(x2, F);
        x3 = x + 0.5*dt*f2;
        f3 = tugboat3d(x3, F);
        x4 = x + dt*f3;
        f4 = tugboat3d(x4, F);
        x  = x + (dt/6)*(f1 + 2*f2 + 2*f3 + f4);
        x(3) = wrapToPiLocal(x(3));
    end

    % Results: trim logs; viol = worst excursion outside the active funnel
    idx = 1:last;
    o.time = time(idx); o.X = X_h(idx); o.Y = Y_h(idx); o.h = h_h(idx);
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

function [un, wn] = durmazNominal(x, ctr, Kv, Ka, rho_tol)
    [rho, alpha] = polarState(x, ctr);
    if rho <= rho_tol
        un = 0; wn = 0;
    else
        un = 2*Kv*rho*cos(alpha);
        wn = Ka*alpha + (un/rho)*sin(alpha);
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

function [M, G, Cn, Dn] = plantMats(nu)
% Reproduces tugboat3d.m's M, C(nu), D(nu) construction exactly.
    m = 10.2; Iz = 0.63994; d = 0.29/2;
    Xu = -5.76909; Xuu = -2.17161; Xud = -0.87818;
    Yv = -3.98659; Yr = -0.0001; Yvv = -3.95131; Yvd = -1.05279; Yrd = -1.92760;
    Nv = -0.0001; Nr = -0.12392; Nrr = -0.33077; Nrd = -0.04531;

    M = [m - Xud, 0,        0;
         0,       m - Yvd, -Yrd;
         0,      -Yrd,      Iz - Nrd];
    B = [1, 1; 0, 0; d, -d];
    G = M \ B;

    u = nu(1); v = nu(2); r = nu(3);
    Cn = [0, 0, -(M(2,2)*v + M(2,3)*r);
          0, 0,   M(1,1)*u;
          M(2,2)*v + M(2,3)*r, -M(1,1)*u, 0];
    Dn = [-Xu - Xuu*abs(u), 0,                 0;
          0,                -Yv - Yvv*abs(v), -Yr;
          0,                -Nv,              -Nr - Nrr*abs(r)];
end

function [b, Lfb, Lf2b, LgLfb] = funnelLie(x, ctr, R)
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
    nu = x(4:6);
    [M, G, Cn, Dn] = plantMats(nu);
    d = -(M \ (Cn*nu + Dn*nu));                 % drift acceleration of nu

    [rho, alpha] = polarState(x, ctr);
    u = nu(1); v = nu(2); r = nu(3);
    wt = -u*sin(alpha) + v*cos(alpha);

    b     = R - rho;
    Lfb   = u*cos(alpha) + v*sin(alpha);
    Lf2b  = cos(alpha)*d(1) + sin(alpha)*d(2) - wt^2/rho - wt*r;
    LgLfb = cos(alpha)*G(1,:) + sin(alpha)*G(2,:);          % 1x2 row, acts on F
end

function F = hocbfForceFilter(x, F_nom, ctr, R, HF, Fmax)
% 2nd-order HOCBF-QP on thrust:
%     min ||F - F_nom||^2
%     s.t. Lf2 b + LgLf b*F + (k1+k2)*Lf b + k1*k2*b >= 0 ,   |F| <= Fmax
% i.e.  -LgLf b * F <= Lf2 b + (k1+k2) Lf b + k1 k2 b.
    [b, Lfb, Lf2b, LgLfb] = funnelLie(x, ctr, R);
    A = -LgLfb;
    c = Lf2b + (HF.k1 + HF.k2)*Lfb + HF.k1*HF.k2*b;

    F = F_nom;
    if A*F_nom > c                              % nominal thrust violates the constraint
        [Fq, ~, exitflag] = quadprog(2*eye(2), -2*F_nom, A, c, [], [], ...
            [-Fmax; -Fmax], [Fmax; Fmax], F_nom, HF.qp_options);
        if exitflag > 0
            F = Fq;
        end
    end
end

function F = pflForce(x, ctr, g, Fmax)
% SwayAwarePFL.m's exact partial-feedback-linearization funnel law.
    u = x(4); v = x(5); r = x(6);
    [rho, alpha] = polarState(x, ctr);
    rho_e = max(rho, g.rho_floor);

    rho_dot   = -u*cos(alpha) - v*sin(alpha);
    alpha_dot = (u*sin(alpha) - v*cos(alpha))/rho_e - r;
    rho_e_dot = (rho > g.rho_floor) * rho_dot;

    uc = g.k_rho * rho * cos(alpha);
    gg = uc*sin(alpha) - v*cos(alpha);
    rc = g.k_alpha*alpha + gg/rho_e;

    uc_dot = g.k_rho * (rho_dot*cos(alpha) - rho*sin(alpha)*alpha_dot);

    eu = u - uc; er = r - rc;

    nu = [u; v; r];
    [M, ~, Cn, Dn] = plantMats(nu);
    m11 = M(1,1); m22 = M(2,2); m23 = M(2,3); m33 = M(3,3);
    fvec = Cn*nu + Dn*nu;
    f1 = fvec(1); f2 = fvec(2); f3 = fvec(3);

    g_dot_known = uc_dot*sin(alpha) + (uc*cos(alpha) + v*sin(alpha))*alpha_dot;

    rc_A = g.k_alpha*alpha_dot + (g_dot_known*rho_e - gg*rho_e_dot)/rho_e^2;
    rc_B = -cos(alpha)/rho_e;
    rc_C = rc_A - rc_B*f2/m22;
    rc_D = -rc_B*m23/m22;

    denom = 1 - rc_D;
    if abs(denom) < 0.1
        denom = 0.1*sign(denom + (denom == 0));
    end
    r_dot_des = (rc_C - g.k_r*er)/denom;

    X  = m11*(uc_dot - g.k_u*eu) + f1;
    Nc = (m33 - m23^2/m22)*r_dot_des - (m23/m22)*f2 + f3;

    d = 0.29/2;
    FL = (X + Nc/d)/2;
    FR = (X - Nc/d)/2;
    F = [min(max(FL, -Fmax), Fmax); min(max(FR, -Fmax), Fmax)];
end

function a = wrapToPiLocal(a)
    a = mod(a + pi, 2*pi) - pi;
end

function s = tf2str(tf)
    if tf, s = 'reached goal'; else, s = 'did not reach goal'; end
end
