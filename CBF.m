%% KINEMATIC CBF / HOCBF SAFETY FILTER ON (u, omega) -> PI LOOPS -> TUGBOAT
%
% Approach 1 (fully kinematic): the Durmaz2024 law produces (u_nom, w_nom).
% A QP filters that pair BEFORE it reaches the low-level PI loops
% (lowLevelControl.m), so the plant and PI are untouched:
%
%   Durmaz2024 (u_nom, w_nom) -> QP filter -> (u_ref, r_ref) -> PI -> tugboat3d
%
% Barrier (single circular funnel, centre c, radius R):
%   h = R - rho,  rho = ||p - c||           (h >= 0 means inside)
%   alpha = phi - psi,  phi = bearing to c, w_t = -u*sin(alpha) + v*cos(alpha)
%
% Exact kinematics with sway (u = surge, v = SWAY, both measured):
%   h_dot  = u*cos(alpha) + v*sin(alpha)
%   h_ddot = u_dot*cos(alpha) + v_dot*sin(alpha) - w_t^2/rho - w_t*r
%
% Constraint A -- CBF on surge, relative degree 1 (u appears in h_dot):
%   u_ref*cos(alpha) + v*sin(alpha) + k1*h >= 0
%   v is READ from the state, not commanded. Without the v term this is
%   u*cos(alpha) >= -k1*h, which Durmaz's law satisfies by construction, so
%   the filter would never act.
%
% Constraint B -- HOCBF on yaw rate, relative degree 2 (r appears in h_ddot):
%   h_ddot + (k1+k2)*h_dot + k1*k2*h >= 0, with the velocity accelerations
%   u_dot, v_dot NEUTRAL (set to 0: this is the kinematic-level model, where
%   velocities are inputs/quasi-static), and measured u, v in w_t and h_dot:
%     -w_t^2/rho - w_t*w_ref + (k1+k2)*h_dot + k1*k2*h >= 0
%   B gives turning authority exactly where A loses it: near the tangent
%   (cos(alpha) ~ 0) surge is useless but sin(alpha) ~ 1. LIMITATION: the
%   w_ref coefficient is -w_t = u*sin(alpha) - v*cos(alpha), which vanishes
%   when the boat is sliding sideways with u ~ 0. Neither A nor B can help
%   there -- only sway ACCELERATION can, which is dynamic (see HOCBF.m).
%
% Both constraints are softened with slacks (heavily penalized) so the QP is
% always feasible; slack use and the "tangent" steps (surge ineffective) are
% logged, because they mark where the kinematic filter is not enough.
%
%   min ||[u_ref;w_ref] - [u_nom;w_nom]||^2 + W*(dA^2 + dB^2)
%
% Scenario is identical to test_durmaz_sideslip.m; each drift angle is run
% twice, Durmaz alone and Durmaz + filter.

clear; clc; close all;

%% ---------------- Simulation setup ----------------
P.dt_sim   = 0.01;                % [s]
P.sim_time = 60.0;                % [s]
P.time     = 0:P.dt_sim:P.sim_time;

%% ---------------- Single funnel ----------------
P.c    = [0; 0];
P.R    = 5.0;                     % funnel radius [m]
P.rho0 = 0.9*P.R;                 % initial distance to center [m]

%% ---------------- Durmaz2024 gains (same as Durmaz2024.m) ----------------
P.Kv = 0.05;
P.Ka = 0.30;
P.rho_arrival_tol = 0.05;

%% ---------------- Low-level control (copied from lowLevelControl.m) ------
% Keep in sync with lowLevelControl.m.
P.beam  = 0.29;
P.d     = P.beam/2;
P.F_max = 500000.0;
P.F_min = -P.F_max;

P.P_yaw = 5.0;   P.I_yaw = 0.02;  P.D_yaw = 0.0;
P.P_speed = 100.0;  P.I_speed = 50.0;  P.D_speed = 0.0;

P.Iyaw_max   = 20000.0;
P.Ispeed_max = 2000*P.F_max;

%% ---------------- Safety-filter parameters (new, tunable) ----------------
% k1 >= U0/h0 (= 4 here) puts every scenario's initial state inside the
% safe set h_dot + k1*h >= 0, which the forward-invariance argument needs.
P.k1 = 5.0;                       % class-K gain on h      [1/s]
P.k2 = 5.0;                       % class-K gain on psi_1  [1/s]
% u_lim is SENSITIVE: at 3 m/s the filter demands full surge and yaw, which
% itself generates large sway (Coriolis m11*u*r) that constraint A cannot
% cancel, and the boat gets trapped orbiting just outside R. At 1 m/s it
% recovers and converges in the whole sweep.
P.u_lim = 1.0;                    % |u_ref| bound [m/s]
P.w_lim = pi/2;                   % |r_ref| bound [rad/s]
P.slack_w = 1e4;                  % slack penalty
P.use_omega_hocbf = true;         % false: constraint A only (surge CBF)
P.tan_tol = 0.1;                  % |cos(alpha)| below this counts as "tangent"
P.leave_tol = 1e-3;               % rho > R + leave_tol counts as leaving (RK4/ZOH discretization) [m]

P.qp_options = optimoptions('quadprog', 'Display', 'off');

%% ---------------- Test scenario ----------------
P.U0     = 2.0;                   % initial speed through water [m/s]
P.alpha0 = pi/2;                  % initial heading error to the center [rad]
beta0_deg = [0 -15 -30 -45 -60 -75 -90];

%% ---------------- Run sweep ----------------
nCase = numel(beta0_deg);
off = cell(1, nCase);
on  = cell(1, nCase);
for i = 1:nCase
    off{i} = runCase(deg2rad(beta0_deg(i)), P, false);
    on{i}  = runCase(deg2rad(beta0_deg(i)), P, true);
end

%% ---------------- Report ----------------
fprintf('\n--- Kinematic CBF filter (k1=%.1f, k2=%.1f, omega-HOCBF %s), R = %.2f m, rho0 = %.2f m ---\n', ...
        P.k1, P.k2, yn(P.use_omega_hocbf), P.R, P.rho0);
fprintf('%10s | %10s %8s | %9s %6s | %8s | %9s | %8s | %8s | %8s\n', ...
        'beta0[deg]', 'maxrho off', 'leaves', 'maxrho on', 'leaves', 'min h on', 'final rho', 'active %', 'tangent', 'slack');
fprintf('%s\n', repmat('-', 1, 117));
for i = 1:nCase
    a = off{i};  b = on{i};
    fprintf('%10.0f | %10.3f %8s | %9.3f %6s | %8.3f | %9.3f | %8.2f | %8d | %8d\n', ...
        beta0_deg(i), a.max_rho, yn(a.left), b.max_rho, yn(b.left), ...
        min(P.R - b.rho), b.rho(end), 100*mean(b.active), sum(b.tangent), sum(b.slack));
end
fprintf('(active %% = steps where the filter changed the command; tangent = steps with |cos(alpha)|<%.2f and\n', P.tan_tol);
fprintf(' surge constraint violated at nominal; slack = steps where the QP needed slack, i.e. constraints unsatisfiable)\n');

%% ---------------- Plots ----------------
colors = parula(nCase + 1);
th = linspace(0, 2*pi, 400);

figure('Name','CBF - trajectories','Color','w');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for s = 1:2
    nexttile; hold on; axis equal; grid on;
    plot(P.c(1) + P.R*cos(th), P.c(2) + P.R*sin(th), 'r-', 'LineWidth', 2);
    plot(P.c(1), P.c(2), 'r+', 'MarkerSize', 10, 'LineWidth', 1.5);
    for i = 1:nCase
        if s == 1, o = off{i}; else, o = on{i}; end
        plot(o.X, o.Y, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
        plot(o.X(1), o.Y(1), 'o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    end
    xlabel('x [m]'); ylabel('y [m]');
    if s == 1, title('Durmaz2024 alone'); else, title('Durmaz2024 + kinematic CBF'); end
end

figure('Name','CBF - time histories','Color','w');
tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on; grid on;
for i = 1:nCase, plot(P.time, off{i}.rho, 'Color', colors(i,:), 'LineWidth', 1.3); end
yline(P.R, 'r--', 'R', 'LineWidth', 1.2);
ylabel('\rho [m]'); title('Distance to center: Durmaz2024 alone');

nexttile; hold on; grid on;
for i = 1:nCase, plot(P.time, on{i}.rho, 'Color', colors(i,:), 'LineWidth', 1.3); end
yline(P.R, 'r--', 'R', 'LineWidth', 1.2);
ylabel('\rho [m]'); title('Distance to center: Durmaz2024 + kinematic CBF');
legend(compose('\\beta_0 = %d\\circ', beta0_deg), 'Location', 'eastoutside');

nexttile; hold on; grid on;
o = on{end};
yyaxis left;  plot(P.time, o.u_nom, '--', P.time, o.u_ref, '-', 'LineWidth', 1.2); ylabel('u [m/s]');
yyaxis right; plot(P.time, o.w_nom, '--', P.time, o.w_ref, '-', 'LineWidth', 1.2); ylabel('r [rad/s]');
xlabel('Time [s]');
title(sprintf('Filter action, \\beta_0 = %d\\circ  (dashed = nominal, solid = filtered)', beta0_deg(end)));

%% ---------------- Functions ----------------
function o = runCase(beta0, P, filter_on)

    N  = numel(P.time);
    dt = P.dt_sim;

    X0   = P.c(1) + P.rho0;
    Y0   = P.c(2);
    phi0 = atan2(P.c(2) - Y0, P.c(1) - X0);
    psi0 = wrapToPiLocal(phi0 - P.alpha0);
    x = [X0; Y0; psi0; P.U0*cos(beta0); P.U0*sin(beta0); 0; 0; 0];

    x_hist   = zeros(8, N);
    rho_hist = zeros(1, N);
    alp_hist = zeros(1, N);
    u_nom_h  = zeros(1, N);  w_nom_h = zeros(1, N);
    u_ref_h  = zeros(1, N);  w_ref_h = zeros(1, N);
    active   = false(1, N);
    tangent  = false(1, N);
    slack    = false(1, N);
    x_hist(:,1) = x;

    int_e_r = 0.0;
    int_e_speed = 0.0;

    for k = 1:N
        psi = x(3);  ub = x(4);  r = x(6);

        % ---- Durmaz2024 nominal law (Eq. 33)
        ex  = x(1) - P.c(1);
        ey  = x(2) - P.c(2);
        rho = hypot(ex, ey);
        alpha = wrapToPiLocal(atan2(-ey, -ex) - psi);

        if rho <= P.rho_arrival_tol
            u_nom = 0;
            w_nom = 0;
        else
            u_nom = 2*P.Kv*rho*cos(alpha);
            w_nom = P.Ka*alpha + (u_nom/rho)*sin(alpha);
        end

        % ---- kinematic safety filter
        u_ref = u_nom;
        w_ref = w_nom;
        if filter_on && rho > P.rho_arrival_tol
            [u_ref, w_ref, inf_] = kinematicFilter(x, u_nom, w_nom, P);
            active(k)  = abs(u_ref - u_nom) > 1e-4 || abs(w_ref - w_nom) > 1e-4;
            tangent(k) = inf_.tangent;
            slack(k)   = inf_.slack;
        end

        rho_hist(k) = rho;
        alp_hist(k) = alpha;
        u_nom_h(k) = u_nom;  w_nom_h(k) = w_nom;
        u_ref_h(k) = u_ref;  w_ref_h(k) = w_ref;
        if k == N, break; end

        % ---- low-level PI loops (as in lowLevelControl.m)
        e_speed = u_ref - ub;
        e_r     = w_ref - r;

        int_e_speed = min(max(int_e_speed + e_speed*dt, -P.Ispeed_max), P.Ispeed_max);
        int_e_r     = min(max(int_e_r     + e_r    *dt, -P.Iyaw_max  ), P.Iyaw_max  );

        tau_X = P.P_speed*e_speed + P.I_speed*int_e_speed - P.D_speed*0;
        tau_N = P.P_yaw  *e_r     + P.I_yaw  *int_e_r     - P.D_yaw  *0;

        FL_cmd = tau_X/2 + tau_N/(2*P.d);
        FR_cmd = tau_X/2 - tau_N/(2*P.d);

        FL_sat = min(max(FL_cmd, P.F_min), P.F_max);
        FR_sat = min(max(FR_cmd, P.F_min), P.F_max);

        if (FL_sat ~= FL_cmd) || (FR_sat ~= FR_cmd)
            tau_X_ach = FL_sat + FR_sat;
            tau_N_ach = P.d*(FL_sat - FR_sat);
            if P.I_speed > 0
                int_e_speed = int_e_speed + (tau_X_ach - tau_X)/P.I_speed;
                int_e_speed = min(max(int_e_speed, -P.Ispeed_max), P.Ispeed_max);
            end
            if P.I_yaw > 0
                int_e_r = int_e_r + (tau_N_ach - tau_N)/P.I_yaw;
                int_e_r = min(max(int_e_r, -P.Iyaw_max), P.Iyaw_max);
            end
        end

        u_cmd = [FL_sat; FR_sat];

        % ---- plant (RK4, ZOH input)
        f1 = tugboat3d(x,               u_cmd);
        f2 = tugboat3d(x + 0.5*dt*f1,   u_cmd);
        f3 = tugboat3d(x + 0.5*dt*f2,   u_cmd);
        f4 = tugboat3d(x +     dt*f3,   u_cmd);
        x  = x + (dt/6)*(f1 + 2*f2 + 2*f3 + f4);

        x_hist(:,k+1) = x;
    end

    U    = hypot(x_hist(4,:), x_hist(5,:));
    beta = atan2(x_hist(5,:), abs(x_hist(4,:)));
    beta(U < 0.05) = NaN;

    o.X = x_hist(1,:);  o.Y = x_hist(2,:);
    o.rho = rho_hist;   o.alpha = alp_hist;  o.beta = beta;
    o.u_nom = u_nom_h;  o.w_nom = w_nom_h;
    o.u_ref = u_ref_h;  o.w_ref = w_ref_h;
    o.active = active;  o.tangent = tangent;  o.slack = slack;

    o.max_rho = max(rho_hist);
    o.left    = any(rho_hist > P.R + P.leave_tol);
    o.t_out   = dt*sum(rho_hist > P.R + P.leave_tol);
end

function [u_ref, w_ref, info] = kinematicFilter(x, u_nom, w_nom, P)
% QP over z = [u_ref; w_ref; dA; dB]:
%   min ||z(1:2) - [u_nom; w_nom]||^2 + W*(dA^2 + dB^2)
%   A:  -ca*u - dA               <= v*sa + k1*h
%   B:   w_t*w - dB              <= -w_t^2/rho + (k1+k2)*h_dot + k1*k2*h

    k1 = P.k1;  k2 = P.k2;

    ex  = x(1) - P.c(1);
    ey  = x(2) - P.c(2);
    rho = hypot(ex, ey);
    al  = atan2(-ey, -ex) - x(3);
    ca  = cos(al);  sa = sin(al);

    u = x(4);  v = x(5);
    w_t  = -u*sa + v*ca;
    h    = P.R - rho;
    hdot = u*ca + v*sa;

    rhsA = v*sa + k1*h;
    rhsB = -w_t^2/rho + (k1+k2)*hdot + k1*k2*h;

    viol_A = (-ca*u_nom > rhsA);
    viol_B = P.use_omega_hocbf && (w_t*w_nom > rhsB);

    info.tangent = viol_A && abs(ca) < P.tan_tol;
    info.slack   = false;

    u_ref = u_nom;
    w_ref = w_nom;
    if ~viol_A && ~viol_B
        return;
    end

    Wd = P.slack_w;
    H  = diag([2, 2, 2*Wd, 2*Wd]);
    f  = [-2*u_nom; -2*w_nom; 0; 0];

    A = [-ca,  0, -1,  0];
    b = rhsA;
    if P.use_omega_hocbf
        A = [A; 0, w_t, 0, -1];
        b = [b; rhsB];
    end

    lb = [-P.u_lim; -P.w_lim; 0;   0];
    ub = [ P.u_lim;  P.w_lim; inf; inf];

    [z, ~, exitflag] = quadprog(H, f, A, b, [], [], lb, ub, [u_nom; w_nom; 0; 0], P.qp_options);

    if exitflag > 0
        u_ref = z(1);
        w_ref = z(2);
        info.slack = any(z(3:4) > 1e-6);
    else
        info.slack = true;      % solver failure: keep the nominal command
    end
end

function s = yn(tf)
    if tf, s = 'YES'; else, s = 'no'; end
end

function a = wrapToPiLocal(a)
    a = mod(a + pi, 2*pi) - pi;
end
