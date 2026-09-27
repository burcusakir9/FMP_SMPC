%% DYNAMIC HOCBF SAFETY FILTER ON THRUSTER FORCES (tugboat3d, relative degree 2)
%
% Approach 2 (dynamic): the Durmaz2024 law and the low-level PI loops
% (lowLevelControl.m) run UNCHANGED and produce a nominal thrust pair
% F_nom = [F_L; F_R]. A QP then filters the FORCES sent to the plant:
%
%   Durmaz2024 (u_nom, w_nom) -> PI -> F_nom -> HOCBF-QP -> F -> tugboat3d
%
% Barrier (single circular funnel, centre c, radius R):
%   h = R - rho,  alpha = phi - psi,  w_t = -u*sin(alpha) + v*cos(alpha)
%   h_dot  = u*cos(alpha) + v*sin(alpha)                       (no force)
%   h_ddot = u_dot*cos(alpha) + v_dot*sin(alpha) - w_t^2/rho - w_t*r
%
% Plant: nu_dot = M^-1*(B*F - C(nu)nu - D(nu)nu). Writing G = M^-1*B and
% d = -M^-1*(C + D)nu, u_dot = G1*F + d1 and v_dot = G2*F + d2, so
%
%   h_ddot = a(x)*F + b(x),
%   a = cos(alpha)*G1 + sin(alpha)*G2                          (1x2)
%   b = cos(alpha)*d1 + sin(alpha)*d2 - w_t^2/rho - w_t*r
%
% a never vanishes (G is invertible and cos, sin are never both 0), unlike
% the kinematic case where surge authority dies at the tangent. Near the
% tangent (cos ~ 0) a ~ G2 = [-0.07, +0.07]: DIFFERENTIAL thrust, which
% accelerates the vessel in sway. That is exactly the authority the
% kinematic filter (CBF.m) lacks when the boat slides sideways at u ~ 0.
%
% HOCBF (Xiao & Belta; linear class-K functions, gains k1, k2):
%   psi_1 = h_dot + k1*h,   psi_2 = psi_1_dot + k2*psi_1 >= 0
%   =>  a*F + b + (k1+k2)*h_dot + k1*k2*h >= 0
%
% QP (single constraint, box limits on the forces):
%   min ||F - F_nom||^2   s.t.  -a*F <= b + (k1+k2)*h_dot + k1*k2*h,
%                               F_min <= F <= F_max
%
% ACTUATOR LAG: tugboat3d.m integrates the 0.25 s thruster lag state but
% applies the COMMANDED force to the dynamics (line 74: F = u(1:2)). So in
% this simulator F_cmd enters h_ddot directly and the constraint above is
% exact for the simulated plant. If the plant is later changed to apply the
% lagged force, h has relative degree 3 in F_cmd and a third-order CBF is
% needed.
%
% Integrators: the PI is deliberately OBLIVIOUS to the filter (only its own
% saturation anti-windup, as in lowLevelControl.m). Two alternatives were
% tried and were worse: unwinding the integrators to the applied force
% leaves a persistent yaw bias (I_yaw = 0.02 takes minutes to wash it out)
% and stalls the boat ~0.35 m from the center; freezing them while the
% filter is active pins the boat on the boundary and it never leaves.
% Results are therefore sensitive to this PI-vs-filter interplay.
%
% Scenario is identical to test_durmaz_sideslip.m; each drift angle is run
% twice, Durmaz + PI alone and Durmaz + PI + HOCBF filter.

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

%% ---------------- HOCBF parameters (new, tunable) ----------------
% k1 >= U0/h0 (= 4 here) puts every scenario's initial state inside the
% safe set h_dot + k1*h >= 0, which the forward-invariance argument needs.
P.k1 = 5.0;                       % class-K gain on h      [1/s]
P.k2 = 5.0;                       % class-K gain on psi_1  [1/s]
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
fprintf('\n--- Dynamic HOCBF filter on thrust (k1=%.1f, k2=%.1f), R = %.2f m, rho0 = %.2f m ---\n', ...
        P.k1, P.k2, P.R, P.rho0);
fprintf('%10s | %10s %8s | %9s %6s | %8s | %9s | %8s | %10s | %8s\n', ...
        'beta0[deg]', 'maxrho off', 'leaves', 'maxrho on', 'leaves', 'min h on', 'final rho', 'active %', 'peak |dF|', 'QP fails');
fprintf('%s\n', repmat('-', 1, 118));
for i = 1:nCase
    a = off{i};  b = on{i};
    fprintf('%10.0f | %10.3f %8s | %9.3f %6s | %8.4f | %9.3f | %8.2f | %10.1f | %8d\n',...
        beta0_deg(i), a.max_rho, yn(a.left), b.max_rho, yn(b.left), ...
        min(P.R - b.rho), b.rho(end), 100*mean(b.active), b.peak_dF, sum(b.qp_fail));
end
fprintf('(active %% = steps where the QP changed the force; peak |dF| in N; QP fails = solver did not converge, nominal kept)\n');

%% ---------------- Plots ----------------
colors = parula(nCase + 1);
th = linspace(0, 2*pi, 400);

figure('Name','HOCBF - trajectories','Color','w');
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
    if s == 1, title('Durmaz2024 + PI alone'); else, title('Durmaz2024 + PI + dynamic HOCBF'); end
end

figure('Name','HOCBF - time histories','Color','w');
tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on; grid on;
for i = 1:nCase, plot(P.time, off{i}.rho, 'Color', colors(i,:), 'LineWidth', 1.3); end
yline(P.R, 'r--', 'R', 'LineWidth', 1.2);
ylabel('\rho [m]'); title('Distance to center: Durmaz2024 + PI alone');

nexttile; hold on; grid on;
for i = 1:nCase, plot(P.time, on{i}.rho, 'Color', colors(i,:), 'LineWidth', 1.3); end
yline(P.R, 'r--', 'R', 'LineWidth', 1.2);
ylabel('\rho [m]'); title('Distance to center: Durmaz2024 + PI + dynamic HOCBF');
legend(compose('\\beta_0 = %d\\circ', beta0_deg), 'Location', 'eastoutside');

nexttile; hold on; grid on;
o = on{end};
plot(P.time, o.F_nom(1,:), '--', P.time, o.F_app(1,:), '-', ...
     P.time, o.F_nom(2,:), '--', P.time, o.F_app(2,:), '-', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('F [N]');
legend('F_L nominal', 'F_L filtered', 'F_R nominal', 'F_R filtered', 'Location', 'eastoutside');
title(sprintf('Filter action, \\beta_0 = %d\\circ', beta0_deg(end)));

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
    F_nom_h  = zeros(2, N);
    F_app_h  = zeros(2, N);
    active   = false(1, N);
    qp_fail  = false(1, N);
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

        rho_hist(k) = rho;
        alp_hist(k) = alpha;
        if k == N, break; end

        % ---- low-level PI loops (as in lowLevelControl.m)
        e_speed = u_nom - ub;
        e_r     = w_nom - r;

        int_e_speed = min(max(int_e_speed + e_speed*dt, -P.Ispeed_max), P.Ispeed_max);
        int_e_r     = min(max(int_e_r     + e_r    *dt, -P.Iyaw_max  ), P.Iyaw_max  );

        tau_X = P.P_speed*e_speed + P.I_speed*int_e_speed - P.D_speed*0;
        tau_N = P.P_yaw  *e_r     + P.I_yaw  *int_e_r     - P.D_yaw  *0;

        FL_cmd = tau_X/2 + tau_N/(2*P.d);
        FR_cmd = tau_X/2 - tau_N/(2*P.d);

        F_nom = [min(max(FL_cmd, P.F_min), P.F_max);
                 min(max(FR_cmd, P.F_min), P.F_max)];

        % ---- PI's own saturation anti-windup (back-calculation)
        if (F_nom(1) ~= FL_cmd) || (F_nom(2) ~= FR_cmd)
            tau_X_ach = F_nom(1) + F_nom(2);
            tau_N_ach = P.d*(F_nom(1) - F_nom(2));
            if P.I_speed > 0
                int_e_speed = int_e_speed + (tau_X_ach - tau_X)/P.I_speed;
                int_e_speed = min(max(int_e_speed, -P.Ispeed_max), P.Ispeed_max);
            end
            if P.I_yaw > 0
                int_e_r = int_e_r + (tau_N_ach - tau_N)/P.I_yaw;
                int_e_r = min(max(int_e_r, -P.Iyaw_max), P.Iyaw_max);
            end
        end

        % ---- dynamic HOCBF filter on the forces
        F = F_nom;
        if filter_on && rho > P.rho_arrival_tol
            [h, hdot, a, b] = funnelBarrier(x, P.c, P.R);
            K = (P.k1 + P.k2)*hdot + P.k1*P.k2*h;

            if a*F_nom + b + K < 0
                [Fq, ~, exitflag] = quadprog(2*eye(2), -2*F_nom, -a, b + K, [], [], ...
                                             [P.F_min; P.F_min], [P.F_max; P.F_max], F_nom, P.qp_options);
                if exitflag > 0
                    F = Fq;
                else
                    qp_fail(k) = true;
                end
            end
            active(k) = norm(F - F_nom) > 1e-3;
        end

        F_nom_h(:,k) = F_nom;
        F_app_h(:,k) = F;


        % ---- plant (RK4, ZOH input)
        f1 = tugboat3d(x,               F);
        f2 = tugboat3d(x + 0.5*dt*f1,   F);
        f3 = tugboat3d(x + 0.5*dt*f2,   F);
        f4 = tugboat3d(x +     dt*f3,   F);
        x  = x + (dt/6)*(f1 + 2*f2 + 2*f3 + f4);

        x_hist(:,k+1) = x;
    end

    U    = hypot(x_hist(4,:), x_hist(5,:));
    beta = atan2(x_hist(5,:), abs(x_hist(4,:)));
    beta(U < 0.05) = NaN;

    o.X = x_hist(1,:);  o.Y = x_hist(2,:);
    o.rho = rho_hist;   o.alpha = alp_hist;  o.beta = beta;
    o.F_nom = F_nom_h;  o.F_app = F_app_h;
    o.active = active;  o.qp_fail = qp_fail;
    o.peak_dF = max(vecnorm(F_app_h - F_nom_h));

    o.max_rho = max(rho_hist);
    o.left    = any(rho_hist > P.R + P.leave_tol);
    o.t_out   = dt*sum(rho_hist > P.R + P.leave_tol);
end

% >>> funnelBarrier begin
function [h, hdot, a, b] = funnelBarrier(x, c, R)
% x = tugboat3d state, c = funnel centre [2x1], R = radius
% Returns h, h_dot and (a, b) with h_ddot = a*F + b, F = [F_L; F_R].
    m=10.2; Iz=0.63994; d=0.145;
    Xu=-5.76909; Xuu=-2.17161; Xud=-0.87818;
    Yv=-3.98659; Yr=-0.0001; Yvv=-3.95131; Yvd=-1.05279; Yrd=-1.92760;
    Nv=-0.0001; Nr=-0.12392; Nrr=-0.33077; Nrd=-0.04531;
    M=[m-Xud 0 0; 0 m-Yvd -Yrd; 0 -Yrd Iz-Nrd];
    B=[1 1; 0 0; d -d];  G=M\B;
    u=x(4); v=x(5); r=x(6); nu=[u;v;r];
    C=[0 0 -(M(2,2)*v+M(2,3)*r); 0 0 M(1,1)*u; M(2,2)*v+M(2,3)*r -M(1,1)*u 0];
    D=[-Xu-Xuu*abs(u) 0 0; 0 -Yv-Yvv*abs(v) -Yr; 0 -Nv -Nr-Nrr*abs(r)];
    dr=-(M\(C*nu+D*nu));
    xr=x(1)-c(1); yr=x(2)-c(2); rho=hypot(xr,yr);
    al=atan2(-yr,-xr)-x(3); ca=cos(al); sa=sin(al);
    w=-u*sa+v*ca;
    h=R-rho;
    hdot=u*ca+v*sa;
    a=ca*G(1,:)+sa*G(2,:);
    b=ca*dr(1)+sa*dr(2)-w^2/rho-w*r;
end
% >>> funnelBarrier end

function s = yn(tf)
    if tf, s = 'YES'; else, s = 'no'; end
end

function a = wrapToPiLocal(a)
    a = mod(a + pi, 2*pi) - pi;
end
