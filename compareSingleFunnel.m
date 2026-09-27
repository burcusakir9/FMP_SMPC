%% COMPARE PFL vs KINEMATIC CBF vs DYNAMIC HOCBF -- SINGLE FUNNEL
%
% One circular funnel, vessel starts near the boundary with its bow
% tangent to the circle and a sweep of initial sideslip (drift) angles
% (same scenario as test_durmaz_sideslip.m). Compares three ways of
% keeping the Durmaz2024 funnel law safe on the real 3-DOF tugboat3d
% plant:
%
%   Durmaz: plain Durmaz2024 law -> PI loops -> F (no safety filter; baseline)
%   PFL   : SwayAwarePFL.m's own exact partial-feedback-linearization
%           funnel-tracking law (no filter, own funnel controller)
%   CBF   : Durmaz2024 (u,w) -> kinematic CBF/HOCBF filter (CBF.m) ->
%           PI loops (lowLevelControl.m) -> F
%   HOCBF : Durmaz2024 (u,w) -> PI loops -> dynamic HOCBF force filter
%           (HOCBF.m) -> F
%
% ---------------------------------------------------------------------
% SYMBOLS USED IN THIS FILE
%   x              whole 8-entry state vector [X Y psi u v r F_L F_R]
%   X, Y           position in the world frame [m]
%   psi            heading angle (where the bow points) [rad]
%   u, v, r        surge (forward), sway (sideways, + to port), yaw rate  (body frame)
%   F_L, F_R       thruster force states
%   nu = [u v r]   body-frame velocities = x(4:6)
%   F = [F_L; F_R] thrust command (the controller output)
%
%   center, R    funnel center and radius
%   ex, ey       vessel position RELATIVE to the center (X-cx, Y-cy)
%   rho          distance from vessel to center = hypot(ex, ey)
%   phi          bearing from vessel to center (world frame)
%   alpha        heading error to the center = phi - psi
%
%   theta0       where the vessel STARTS around the center (position angle)
%   rho0, alpha0 initial distance / initial heading error to the center
%   U0, beta0    initial speed magnitude and drift angle:  u = U0 cos(beta0),
%                v = U0 sin(beta0)   (U0 is a speed, not a state)
%
%   u_nominal, w_nominal   Durmaz2024 surge / YAW-RATE references
%   u_ref, w_ref           the same after the CBF filter
%   (w_* is a reference for the yaw rate r; r itself is the measured one)
% ---------------------------------------------------------------------

clear; clc; close all;

%% Scenario
R      = 5.0;              % funnel radius [m]
center    = [0; 0];           % funnel center
rho0   = 0.9*R;            % initial distance to center [m]
theta0 = 0;                % where the vessel starts around the center [rad] (position angle, NOT the heading)
alpha0 = pi/2;             % initial heading error to the center [rad] (pi/2: bow tangent to the circle)
U0     = 2.0;              % initial speed magnitude sqrt(u^2+v^2) [m/s]
betas_deg = [0 -15 -30 -45 -60 -75 -90];   % initial drift angles beta0 [deg] (angle between velocity and bow; < 0: sway away from center)

dt   = 0.01;
Tsim = 60;
Fmax = 100;                % thruster force limit [N]

%% Durmaz20x24 nominal law (CBF & HOCBF branches)
Kv = 0.05; 
Ka = 0.30; 
rho_tol = 0.05;

%% Low-level PI (lowLevelControl.m gains)
PI.P_speed = 100; 
PI.I_speed = 50;
PI.P_yaw   = 5;   
PI.I_yaw   = 0.02;
PI.Ispeed_max = 2000*Fmax; 
PI.Iyaw_max = 20000;

%% Kinematic CBF/HOCBF filter on (u,w) refs (CBF.m)
KF.k1 = 5; 
KF.k2 = 5;
KF.u_lim = 1.0; 
KF.w_lim = pi/2;
KF.use_omega_hocbf = true;
KF.slack_w = 1e4;
KF.qp_options = optimoptions('quadprog', 'Display', 'off');

%% Dynamic HOCBF force filter (HOCBF.m)
HF.k1 = 5;
HF.k2 = 5;
HF.qp_options = optimoptions('quadprog', 'Display', 'off');

%% PFL gains (SwayAwarePFL.m defaults)
PFLg.k_rho = 0.10; 
PFLg.k_alpha = 0.30;
PFLg.k_u   = 1.0;  
PFLg.k_r     = 1.0;
PFLg.rho_floor = 1.0;

%% Run the drift-angle sweep
controllers = {'Durmaz', 'PFL', 'CBF', 'HOCBF'};
nB = numel(betas_deg);
res = cell(numel(controllers), nB);
for c = 1:numel(controllers)
    for b = 1:nB
        res{c,b} = runSingleFunnel(controllers{c}, deg2rad(betas_deg(b)), R, center, ...
            rho0, theta0, alpha0, U0, dt, Tsim, Fmax, Kv, Ka, rho_tol, PI, KF, HF, PFLg);
    end
end

%% Report
fprintf('\n--- Single funnel, R = %.1f m, Fmax = %.0f N ---\n', R, Fmax);
fprintf('%10s', 'beta0[deg]'); fprintf('%9d', betas_deg); fprintf('\n');
for c = 1:numel(controllers)
    fprintf('%10s', controllers{c});
    for b = 1:nB, fprintf('%9.3f', res{c,b}.viol); end
    fprintf('\n');
end
fprintf('(excursion outside the funnel radius, in meters; 0 = stayed inside)\n');

%% Plot 1: trajectories
th = linspace(0, 2*pi, 400);
colors = parula(nB+1);
figure('Name','Single funnel - trajectories','Color','w');
tiledlayout(1,4,'Padding','compact','TileSpacing','compact');
for c = 1:numel(controllers)
    nexttile; hold on; axis equal; grid on;
    fill(center(1)+R*cos(th), center(2)+R*sin(th), [1 0.93 0.85], 'EdgeColor',[0.9 0.4 0],'LineWidth',1.4);
    for b = 1:nB
        plot(res{c,b}.X, res{c,b}.Y, 'Color', colors(b,:), 'LineWidth', 1.3);
        plot(res{c,b}.X(1), res{c,b}.Y(1), '.', 'Color', colors(b,:), 'MarkerSize', 10);
    end
    plot(center(1), center(2), 'k+');
    title(controllers{c}); xlabel('x [m]'); if c == 1, ylabel('y [m]'); end
    axis([-8 8 -8 8]);
end
cb = colorbar; cb.Layout.Tile = 'east'; colormap(colors(1:end-1,:));
clim([min(betas_deg) max(betas_deg)]); cb.Label.String = '\beta_0 [deg]';

%% Plot 2: worst case + excursion summary
figure('Name','Single funnel - comparison','Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
for c = 1:numel(controllers)
    plot(res{c,end}.time, res{c,end}.rho, 'LineWidth', 1.4);
end
yline(R, 'k--', 'R');
xlabel('t [s]'); ylabel('\rho [m]');
title(sprintf('Distance to center, \\beta_0 = %d^\\circ', betas_deg(end)));
legend(controllers, 'Location', 'best');

nexttile; hold on; grid on;
for c = 1:numel(controllers)
    plot(-betas_deg, cellfun(@(results) results.viol, res(c,:)), 'o-', 'LineWidth', 1.4);
end
xlabel('|\beta_0| [deg]'); ylabel('max excursion outside funnel [m]');
title('Excursion vs initial drift angle');
legend(controllers, 'Location', 'best');

%% ======================================================  Local functions

% Run single funnel simulation
function results = runSingleFunnel(controller, beta0, R, center, rho0, theta0, alpha0, U0, dt, Tsim, ...
        Fmax, Kv, Ka, rho_tol, PI, KF, HF, PFLg)

    N = round(Tsim/dt) + 1;
    time = (0:N-1)*dt;

    % Initial states
    % position: polar coordinates (rho0, theta0) around the funnel center
    X0 = center(1) + rho0*cos(theta0); 
    Y0 = center(2) + rho0*sin(theta0);
    % bearing from the vessel to the center (= theta0 + pi)
    phi0 = atan2(center(2) - Y0, center(1) - X0);
    % heading psi0 chosen so that the heading error phi - psi equals alpha0
    psi0 = wrapToPiLocal(phi0 - alpha0);

    % Initial state array [X Y psi u v r F_L F_R] (pose, body velocities: surge/sway/yaw rate, thruster force states)
    x = [X0; Y0; psi0; U0*cos(beta0); U0*sin(beta0); 0; 0; 0]; 

    % Logs
    X_h = zeros(1,N); 
    Y_h = zeros(1,N); 
    rho_h = zeros(1,N);
    int_speed = 0; 
    int_yaw = 0;

    % Time loop
    for k = 1:N
        % Distance to funnel center
        [rho, ~] = polarState(x, center);
        X_h(k) = x(1); 
        Y_h(k) = x(2); 
        rho_h(k) = rho;

        if k == N, break; end

        % Compute thrust F = [F_L; F_R] according to the controller
        switch controller
            case 'PFL'
                F = pflForce(x, center, PFLg, Fmax);

            otherwise % 'CBF' or 'HOCBF': Durmaz2024 nominal law + PI
                % Durmaz2024 nominal surge / yaw-rate references
                [u_nominal, w_nominal] = durmazNominal(x, center, Kv, Ka, rho_tol);
                u_ref = u_nominal; 
                w_ref = w_nominal;
                % CBF: correct the references before they reach the PI loops
                if strcmp(controller, 'CBF') && rho > rho_tol
                    [u_ref, w_ref] = kinFilter(x, u_nominal, w_nominal, center, R, KF);
                end

                % PI controller
                % errors: reference minus measured surge x(4)=u and yaw rate x(6)=r
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
                if strcmp(controller, 'HOCBF') && rho > rho_tol
                    F = hocbfForceFilter(x, F, center, R, HF, Fmax);
                end
        end

        % Plant step (RK4, thrust held constant over dt)
        f1 = tugboat3d(x,             F);
        f2 = tugboat3d(x + 0.5*dt*f1, F);
        f3 = tugboat3d(x + 0.5*dt*f2, F);
        f4 = tugboat3d(x +     dt*f3, F);
        x  = x + (dt/6)*(f1 + 2*f2 + 2*f3 + f4);
        x(3) = wrapToPiLocal(x(3));
    end

    % Results: max excursion outside the funnel radius
    results.time = time; results.X = X_h; results.Y = Y_h; results.rho = rho_h;
    results.max_rho = max(rho_h);
    results.viol = max(0, results.max_rho - R);
end

% Convert cartesian states to polar states
function [rho, alpha] = polarState(x, center)
% Vessel state relative to the funnel center.
    % position relative to the center
    ex = x(1) - center(1);
    ey = x(2) - center(2);
    % distance to the center
    rho = hypot(ex, ey);
    % bearing from the vessel to the center (direction of -[ex; ey])
    phi = atan2(-ey, -ex);
    % heading error: angle between the bow and the direction to the center
    alpha = wrapToPiLocal(phi - x(3));
end

% Durmaz2024 kinematic controller for circular funnels
function [u_nominal, w_nominal] = durmazNominal(x, center, Kv, Ka, rho_tol)
    [rho, alpha] = polarState(x, center);
    if rho <= rho_tol
        u_nominal = 0; w_nominal = 0;
    else
        u_nominal = 2*Kv*rho*cos(alpha);
        w_nominal = Ka*alpha + (u_nominal/rho)*sin(alpha);
    end
end

% Kinematic CBF or HOCBF fiter on the references (u_ref, w_ref).
function [u_ref, w_ref] = kinFilter(x, u_nom, w_nom, center, R, KF)

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

    [rho, alpha] = polarState(x, center);
    u = x(4);                                   % measured surge
    v = x(5);                                   % measured sway (drift term, not an input)
    % the INPUT nu = [u_ref; w_ref] here is (surge ref, yaw-rate ref),
    % not the measured velocities
    wt = -u*sin(alpha) + v*cos(alpha);          % tangential velocity (component of motion perpendicular to the bearing line)

    b = R - rho; % abraxas: or R^2 - rho^2
    

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
    % nominal input already safe
    if ~viol1 && ~viol2
        return;  
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

% Calculate plant model matrices
function [M, G, Cn, Dn] = plantMats(nu)
% Reproduces tugboat3d.m's M, C(nu), D(nu) construction exactly.
    m = 10.2; 
    Iz = 0.63994; 
    d = 0.29/2;
    Xu = -5.76909; 
    Xuu = -2.17161; 
    Xud = -0.87818;
    Yv = -3.98659; 
    Yr = -0.0001; 
    Yvv = -3.95131; 
    Yvd = -1.05279; 
    Yrd = -1.92760;
    Nv = -0.0001; 
    Nr = -0.12392; 
    Nrr = -0.33077; 
    Nrd = -0.04531;

    M = [m - Xud, 0,        0;
         0,       m - Yvd, -Yrd;
         0,      -Yrd,      Iz - Nrd];

    B = [1, 1; 0, 0; d, -d];

    G = M \ B;

    u = nu(1); 
    v = nu(2); 
    r = nu(3);

    Cn = [0, 0, -(M(2,2)*v + M(2,3)*r);
          0, 0,   M(1,1)*u;
          M(2,2)*v + M(2,3)*r, -M(1,1)*u, 0];

    Dn = [-Xu - Xuu*abs(u), 0,                 0;
          0,                -Yv - Yvv*abs(v), -Yr;
          0,                -Nv,              -Nr - Nrr*abs(r)];
end

% Dynamic HOCBF filter
function F = hocbfForceFilter(x, F_nom, center, R, HF, Fmax)
% 2nd-order HOCBF-QP on thrust:
%     min ||F - F_nom||^2
%     s.t. Lf2 b + LgLf b*F + (k1+k2)*Lf b + k1*k2*b >= 0 ,   |F| <= Fmax
% i.e.  -LgLf b * F <= Lf2 b + (k1+k2) Lf b + k1 k2 b.
    [b, Lfb, Lf2b, LgLfb] = funnelLie(x, center, R);
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


function [b, Lfb, Lf2b, LgLfb] = funnelLie(x, center, R)
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
    nu = x(4:6);                                % measured body velocities [u; v; r]
    [M, G, Cn, Dn] = plantMats(nu);
    d = -(M \ (Cn*nu + Dn*nu));                 % drift acceleration of nu

    [rho, alpha] = polarState(x, center);
    u = nu(1); v = nu(2); r = nu(3);
    wt = -u*sin(alpha) + v*cos(alpha);

    b     = R - rho;
    Lfb   = u*cos(alpha) + v*sin(alpha);
    Lf2b  = cos(alpha)*d(1) + sin(alpha)*d(2) - wt^2/rho - wt*r;
    LgLfb = cos(alpha)*G(1,:) + sin(alpha)*G(2,:);          % 1x2 row, acts on F
end

% Sway aware partial deefback limnearization controller
function F = pflForce(x, center, g, Fmax)
% SwayAwarePFL.m's exact partial-feedback-linearization funnel law.
    u = x(4); v = x(5); r = x(6);               % measured surge, sway, yaw rate
    [rho, alpha] = polarState(x, center);
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

% wrap angle between [-pi pi]
function a = wrapToPiLocal(a)
    a = mod(a + pi, 2*pi) - pi;
end
