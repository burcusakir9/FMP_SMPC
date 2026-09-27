%% TEST: DURMAZ2024 FUNNEL LAW ON THE 3-DOF TUGBOAT UNDER LARGE SIDESLIP
%
% Durmaz2024's proof that rho never grows (rho_dot = -2*Kv*rho*cos^2(alpha)
% <= 0, Proposition 1) assumes the KINEMATIC UNICYCLE: velocity always
% points along the heading. The tugboat has sway, so the true range rate is
%
%   rho_dot = -U*cos(alpha - beta),   U = hypot(u,v),  beta = drift angle,
%
% which is only equal to Durmaz's expression when beta = 0. This script
% checks whether a large sideslip lets the vessel leave the funnel.
%
% Architecture (single circular funnel, no RSC tree needed):
%   Durmaz2024 law (v_cmd, omega_cmd)
%     -> lowLevelControl PI loops: u_ref = v_cmd, r_ref = omega_cmd
%     -> thrust allocation -> tugboat3d (RK4, ZOH)
%
% Test scenario (identical for every run except beta0):
%   - vessel starts at rho0 = 0.9*R, inside the funnel
%   - heading tangent to the funnel circle (alpha0 = +90 deg), so Durmaz
%     commands v_cmd ~ 0 and pure yaw toward the center
%   - initial body velocity U0*[cos(beta0); sin(beta0)]. Body y points to
%     port and the center lies to port (alpha0 > 0), so beta0 < 0 is sway
%     directed AWAY from the center. beta0 = 0 is the kinematic baseline.
%
% A run "violates" if rho exceeds R (leaves the funnel) or exceeds rho0
% (contradicts Proposition 1's monotone-rho guarantee, even if still inside).

clear; clc; close all;

%% ---------------- Simulation setup ----------------
P.dt_sim   = 0.01;                % [s]
P.sim_time = 60.0;                % [s]
P.time     = 0:P.dt_sim:P.sim_time;

%% ---------------- Single funnel ----------------
P.c    = [0; 0];                  % funnel center
P.R    = 5.0;                     % funnel radius [m]
P.rho0 = 0.9*P.R;                 % initial distance to center [m]

%% ---------------- Durmaz2024 gains (same as Durmaz2024.m) ----------------
P.Kv = 0.05;
P.Ka = 0.30;
P.rho_arrival_tol = 0.05;         % avoid v/rho singularity at the center

%% ---------------- Low-level control (copied from lowLevelControl.m) ------
% Keep in sync with lowLevelControl.m.
P.beam  = 0.29;
P.d     = P.beam/2;
P.F_max = 500000.0;
P.F_min = -P.F_max;

P.P_yaw = 30.0;   P.I_yaw = 20.0;  P.D_yaw = 0.0;
P.P_speed = 100.0;  P.I_speed = 50.0;  P.D_speed = 0.0;

P.Iyaw_max   = 20000.0;
P.Ispeed_max = 2000*P.F_max;

%% ---------------- Test scenario ----------------
P.U0     = 2.0;                   % initial speed through water [m/s]
P.alpha0 = pi/2;                  % initial heading error to the center [rad]
beta0_deg = [0 -15 -30 -45 -60 -75 -90];   % initial drift angles [deg]

%% ---------------- Run sweep ----------------
nCase = numel(beta0_deg);
res   = cell(1, nCase);
for i = 1:nCase
    res{i} = runCase(deg2rad(beta0_deg(i)), P);
end

%% ---------------- Report ----------------
fprintf('\n--- Durmaz2024 funnel test under sideslip (R = %.2f m, rho0 = %.2f m, U0 = %.1f m/s) ---\n', ...
        P.R, P.rho0, P.U0);
fprintf('%10s | %8s | %11s | %10s | %10s | %14s | %9s\n', ...
        'beta0[deg]', 'max rho', 'rho > rho0', 'leaves R?', 't_out [s]', 'beta@exit[deg]', 'final rho');
fprintf('%s\n', repmat('-', 1, 92));
for i = 1:nCase
    o = res{i};
    fprintf('%10.0f | %8.3f | %11s | %10s | %10.2f | %14.1f | %9.3f\n', ...
        beta0_deg(i), o.max_rho, yn(o.grew), yn(o.left), o.t_out, ...
        rad2deg(o.beta_exit), o.rho(end));
end

%% ---------------- Plots ----------------
colors = parula(nCase + 1);
th = linspace(0, 2*pi, 400);

figure('Name','Durmaz2024 sideslip test - trajectories','Color','w');
hold on; axis equal; grid on;
plot(P.c(1) + P.R*cos(th), P.c(2) + P.R*sin(th), 'r-', 'LineWidth', 2);
plot(P.c(1), P.c(2), 'r+', 'MarkerSize', 10, 'LineWidth', 1.5);
h = gobjects(1, nCase);
for i = 1:nCase
    h(i) = plot(res{i}.X, res{i}.Y, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    plot(res{i}.X(1), res{i}.Y(1), 'o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
end
xlabel('x [m]'); ylabel('y [m]');
title('Single funnel: trajectories for different initial drift angles');
legend(h, compose('\\beta_0 = %d\\circ', beta0_deg), 'Location', 'bestoutside');

figure('Name','Durmaz2024 sideslip test - time histories','Color','w');
tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on; grid on;
for i = 1:nCase
    plot(P.time, res{i}.rho, 'Color', colors(i,:), 'LineWidth', 1.3);
end
yline(P.R,    'r--', 'R',    'LineWidth', 1.2);
yline(P.rho0, 'k:',  '\rho_0', 'LineWidth', 1.0);
ylabel('\rho [m]'); title('Distance to funnel center');

nexttile; hold on; grid on;
for i = 1:nCase
    plot(P.time, rad2deg(res{i}.beta), 'Color', colors(i,:), 'LineWidth', 1.3);
end
ylabel('\beta [deg]'); title('Drift angle \beta = atan2(v, |u|)');

nexttile; hold on; grid on;
for i = 1:nCase
    plot(P.time, rad2deg(res{i}.alpha), 'Color', colors(i,:), 'LineWidth', 1.3);
end
ylabel('\alpha [deg]'); xlabel('Time [s]'); title('Heading error to funnel center');

%% ---------------- Functions ----------------
function o = runCase(beta0, P)

    N = numel(P.time);
    dt = P.dt_sim;

    % Initial state on the +x side of the center, heading tangent.
    X0   = P.c(1) + P.rho0;
    Y0   = P.c(2);
    phi0 = atan2(P.c(2) - Y0, P.c(1) - X0);          % bearing to center (= pi)
    psi0 = wrapToPiLocal(phi0 - P.alpha0);
    x = [X0; Y0; psi0; P.U0*cos(beta0); P.U0*sin(beta0); 0; 0; 0];

    x_hist   = zeros(8, N);
    rho_hist = zeros(1, N);
    alp_hist = zeros(1, N);
    x_hist(:,1) = x;

    int_e_r = 0.0;
    int_e_speed = 0.0;

    for k = 1:N
        psi = x(3);  ub = x(4);  vb = x(5);  r = x(6);

        % ---- Durmaz2024 outer law (Eq. 33), single active funnel
        ex  = x(1) - P.c(1);
        ey  = x(2) - P.c(2);
        rho = hypot(ex, ey);
        phi = atan2(-ey, -ex);
        alpha = wrapToPiLocal(phi - psi);

        if rho <= P.rho_arrival_tol
            v_cmd = 0;
            w_cmd = 0;
        else
            v_cmd = 2*P.Kv*rho*cos(alpha);
            w_cmd = P.Ka*alpha + (v_cmd/rho)*sin(alpha);
        end

        rho_hist(k) = rho;
        alp_hist(k) = alpha;
        if k == N, break; end

        % ---- low-level PI loops (as in lowLevelControl.m)
        e_speed = v_cmd - ub;
        e_r     = w_cmd - r;

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

    % Drift angle, undefined when (nearly) at rest
    U    = hypot(x_hist(4,:), x_hist(5,:));
    beta = atan2(x_hist(5,:), abs(x_hist(4,:)));
    beta(U < 0.05) = NaN;

    o.X     = x_hist(1,:);
    o.Y     = x_hist(2,:);
    o.rho   = rho_hist;
    o.alpha = alp_hist;
    o.beta  = beta;

    o.max_rho      = max(rho_hist);
    o.grew         = o.max_rho > P.rho0 + 1e-6;
    o.left         = any(rho_hist > P.R);
    o.t_out        = dt*sum(rho_hist > P.R);

    % Drift angle at the first funnel exit (NaN if it never leaves)
    k_exit = find(rho_hist > P.R, 1);
    if isempty(k_exit)
        o.beta_exit = NaN;
    else
        o.beta_exit = beta(k_exit);
    end
end

function s = yn(tf)
    if tf, s = 'YES'; else, s = 'no'; end
end

function a = wrapToPiLocal(a)
    a = mod(a + pi, 2*pi) - pi;
end
