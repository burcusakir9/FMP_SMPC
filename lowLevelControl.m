clear; clc; close all;

%% ---------------- Simulation setup ----------------
sim_time = 50.0;              % [s]
dt_sim   = 0.01;               % [s]  (integration step)
time     = 0:dt_sim:sim_time;
N        = numel(time);

%% ---------------- Plant geometry / limits ----------------
% must match tugboat3d.m
beam  = 0.29;                  % [m]
d     = beam/2;                % [m] thruster moment arm
F_max = 500000.0;                  % [N] per-thruster saturation  <-- tune/limit
F_min = -F_max;                % set to 0 for forward-only thrusters

%% ---------------- Controller gains (TUNE HERE) ----------------
% Yaw-rate loop (PI on yaw-rate error; D would need r_dot, kept 0)
P_yaw = 30.0;
I_yaw = 20.0;
D_yaw = 0.0;

% Speed loop (PI on surge speed error)
P_speed = 100.0;
I_speed = 50.0;
D_speed = 0.0;                 % on -du/dt via measured accel-free form (kept 0)

% Integrator clamps (anti-windup)
Iyaw_max   = 20000.0;             % [N*m]
Ispeed_max = 2000*F_max;          % [N]

%% ---------------- References ----------------
speed_ref = 3.0;               % [m/s] surge speed reference
u_ref     = speed_ref * ones(1,N);

r_ref_deg          = zeros(1,N);
r_ref_deg(time >= 20) = 20;                    % [deg/s] yaw-rate reference
r_ref              = deg2rad(r_ref_deg);       % NOTE: deg/s -> rad/s


%% ---------------- Allocate histories ----------------
x_hist    = zeros(8, N);       % plant state
u_hist    = zeros(2, N);       % [F_L_cmd; F_R_cmd]
tau_hist  = zeros(2, N);       % [tau_X; tau_N] requested (pre-saturation)
e_hist    = zeros(2, N);       % [e_speed; e_r]
U_hist    = zeros(1, N);       % total speed sqrt(u^2+v^2)

%% ---------------- Initial condition ----------------
x = zeros(8,1);                % at origin, at rest, thrusters off
% x(3) = deg2rad(0);           % initial heading, if you want an offset

int_e_r     = 0.0;
int_e_speed = 0.0;

x_hist(:,1) = x;
U_hist(1)   = hypot(x(4), x(5));

%% ---------------- Main loop ----------------
for k = 1:N-1

    % ---- measurements
    ub  = x(4);
    vb  = x(5);
    r   = x(6);

    % ---- errors
    e_speed = u_ref(k) - ub;
    e_r     = r_ref(k) - r;

    % ---- integrators (trapezoid-free, simple forward; clamped)
    int_e_speed = min(max(int_e_speed + e_speed*dt_sim, -Ispeed_max), Ispeed_max);
    int_e_r     = min(max(int_e_r     + e_r    *dt_sim, -Iyaw_max  ), Iyaw_max  );

    % ---- control laws
    tau_X = P_speed*e_speed + I_speed*int_e_speed - D_speed*0;   % [N]
    tau_N = P_yaw  *e_r     + I_yaw  *int_e_r     - D_yaw  *0;   % [N*m]

    % ---- thrust allocation (inverse of B_prop rows 1 and 3)
    FL_cmd = tau_X/2 + tau_N/(2*d);
    FR_cmd = tau_X/2 - tau_N/(2*d);

    % ---- saturation
    FL_sat = min(max(FL_cmd, F_min), F_max);
    FR_sat = min(max(FR_cmd, F_min), F_max);

    % ---- back-calculation anti-windup: unwind integrators if clipped
    if (FL_sat ~= FL_cmd) || (FR_sat ~= FR_cmd)
        tau_X_ach = FL_sat + FR_sat;
        tau_N_ach = d*(FL_sat - FR_sat);
        if I_speed > 0
            int_e_speed = int_e_speed + (tau_X_ach - tau_X)/I_speed;
            int_e_speed = min(max(int_e_speed, -Ispeed_max), Ispeed_max);
        end
        if I_yaw > 0
            int_e_r = int_e_r + (tau_N_ach - tau_N)/I_yaw;
            int_e_r = min(max(int_e_r, -Iyaw_max), Iyaw_max);
        end
    end

    u_cmd = [FL_sat; FR_sat];

    % ---- log (control applied over [t_k, t_k+1))
    u_hist(:,k)   = u_cmd;
    tau_hist(:,k) = [tau_X; tau_N];
    e_hist(:,k)   = [e_speed; e_r];

    % ---- integrate plant (RK4, ZOH input)
    f1 = tugboat3d(x,                 u_cmd);
    f2 = tugboat3d(x + 0.5*dt_sim*f1, u_cmd);
    f3 = tugboat3d(x + 0.5*dt_sim*f2, u_cmd);
    f4 = tugboat3d(x +     dt_sim*f3, u_cmd);
    x  = x + (dt_sim/6)*(f1 + 2*f2 + 2*f3 + f4);

    x_hist(:,k+1) = x;
    U_hist(k+1)   = hypot(x(4), x(5));
end

% hold last sample so plots don't drop to zero
u_hist(:,N)   = u_hist(:,N-1);
tau_hist(:,N) = tau_hist(:,N-1);
e_hist(:,N)   = [u_ref(N) - x_hist(4,N); r_ref(N) - x_hist(6,N)];

%% ---------------- Unpack ----------------
X   = x_hist(1,:);  Y   = x_hist(2,:);  psi = x_hist(3,:);
ub  = x_hist(4,:);  vb  = x_hist(5,:);  r   = x_hist(6,:);
FL  = x_hist(7,:);  FR  = x_hist(8,:);

%% ---------------- Performance metrics ----------------
fprintf('\n--- Performance ---\n');
fprintf('Final surge speed      : %.3f m/s (ref %.3f)\n', ub(end), speed_ref);
fprintf('Final yaw rate         : %.2f deg/s (ref %.2f)\n', rad2deg(r(end)), r_ref_deg(end));
fprintf('Steady-state speed err : %.4f m/s\n', e_hist(1,end));
fprintf('Steady-state yaw-rate err : %.4f deg/s\n', rad2deg(e_hist(2,end)));
fprintf('RMS yaw-rate error     : %.4f deg/s\n', rad2deg(sqrt(mean(e_hist(2,:).^2))));
fprintf('RMS speed error        : %.4f m/s\n', sqrt(mean(e_hist(1,:).^2)));
fprintf('Max |sway|             : %.4f m/s\n', max(abs(vb)));
fprintf('Thrust saturation      : %.1f %% of samples\n', ...
        100*mean(abs(u_hist(1,:)) >= F_max-1e-9 | abs(u_hist(2,:)) >= F_max-1e-9));

%% ---------------- Plots ----------------
figure('Name','Tugboat PID - tracking','Color','w');

subplot(3,2,1);
plot(Y, X, 'LineWidth', 1.3); hold on;
plot(Y(1), X(1), 'go', 'MarkerFaceColor','g');
plot(Y(end), X(end), 'rs', 'MarkerFaceColor','r');
xlabel('East, Y [m]'); ylabel('North, X [m]');
title('Trajectory'); grid on; axis equal;

subplot(3,2,2);
plot(time, rad2deg(r), 'LineWidth', 1.3); hold on;
plot(time, r_ref_deg, 'k--', 'LineWidth', 1.0);
xlabel('t [s]'); ylabel('r [deg/s]');
title('Yaw rate'); legend('r','r_{ref}','Location','best'); grid on;

subplot(3,2,3);
plot(time, ub, 'LineWidth', 1.3); hold on;
plot(time, U_hist, 'LineWidth', 1.0);
plot(time, u_ref, 'k--', 'LineWidth', 1.0);
xlabel('t [s]'); ylabel('speed [m/s]');
title('Surge speed'); legend('u','U = |\nu_{1:2}|','u_{ref}','Location','best'); grid on;

subplot(3,2,4);
yyaxis left;  plot(time, rad2deg(e_hist(2,:)), 'LineWidth', 1.2); ylabel('e_r [deg/s]');
yyaxis right; plot(time, e_hist(1,:), 'LineWidth', 1.2);           ylabel('e_u [m/s]');
xlabel('t [s]'); title('Tracking errors'); grid on;

subplot(3,2,5);
plot(time, u_hist(1,:), '--', 'LineWidth', 1.0); hold on;
plot(time, u_hist(2,:), '--', 'LineWidth', 1.0);
plot(time, FL, 'LineWidth', 1.3);
plot(time, FR, 'LineWidth', 1.3);
yline( F_max, 'k:'); yline(F_min, 'k:');
xlabel('t [s]'); ylabel('F [N]');
title('Thruster commands vs. actual (1st-order actuator)');
legend('F_{L,cmd}','F_{R,cmd}','F_L','F_R','Location','best'); grid on;

subplot(3,2,6);
yyaxis left;  plot(time, vb, 'LineWidth', 1.2);            ylabel('v [m/s]');
yyaxis right; plot(time, rad2deg(psi), 'LineWidth', 1.2);  ylabel('\psi [deg]');
xlabel('t [s]'); title('Sway and heading'); grid on;

sgtitle(sprintf('Tugboat PI (yaw rate + speed):  P_{r}=%.3g, I_{r}=%.3g, D_{r}=%.3g  |  P_{u}=%.3g, I_{u}=%.3g', ...
        P_yaw, I_yaw, D_yaw, P_speed, I_speed));

%% ---------------- Helper ----------------
function a = wrapToPiLocal(a)
    a = mod(a + pi, 2*pi) - pi;
end