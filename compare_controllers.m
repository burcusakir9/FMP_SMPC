%% CONTROLLER COMPARISON: Ege & Ankarali (2019) vs Durmaz, Ozdemir & Ankarali (2024)
%
% Both controllers steer a first-order unicycle through the SAME funnel
% chain, using the same funnel-selection scaffold. Only the inner control
% law differs, so building one funnel tree and running both scripts on
% it isolates exactly that difference:
%
%   Ege2019    : v = Kv*rho*cos(phi),      omega = Kphi*phi
%   Durmaz2024 : v = 2*Kv*rho*cos(alpha),  omega = Ka*alpha + (v/rho)*sin(alpha)
%
% Durmaz2024's extra (v/rho)*sin(alpha) term is a feedback-linearizing
% correction that provably keeps rho non-increasing (Proposition 1 in
% the 2024 paper), i.e. the vehicle can't leave the active funnel --
% a guarantee Ege2019's simpler policy does not have.
%
% Note: each script uses its own hardcoded gains (not velocity-matched
% between the two controllers), so this compares the two methods as
% they ship, not gain-for-gain.

close all; clear; clc;

set(0, 'DefaultFigureVisible', 'off'); % suppress each script's own figures

%% BUILD SHARED FUNNEL TREE

fprintf('Building shared funnel tree (RSC.m)...\n');
run('RSC.m');

if isempty(pathIds)
    error('RSC produced no path -- nothing to compare controllers on.');
end

fprintf('Funnel tree built: %d nodes, %d-funnel path.\n\n', numel(nodes), numel(pathIds));

%% RUN EGE 2019

fprintf('Running Ege & Ankarali (2019)...\n');
run('Ege2019.m');
results.ege = packResults(time, state_hist, v_hist, omega_hist, q_goal, goal_tol);
fprintf('  %s\n\n', summaryLine(results.ege));

close all;

%% RUN DURMAZ 2024

fprintf('Running Durmaz, Ozdemir & Ankarali (2024)...\n');
run('Durmaz2024.m');
results.durmaz = packResults(time, state_hist, v_hist, omega_hist, q_goal, goal_tol);
fprintf('  %s\n\n', summaryLine(results.durmaz));

close all;

%% METRICS TABLE

fprintf('============================================================\n');
fprintf('                 CONTROLLER COMPARISON\n');
fprintf('============================================================\n');
fprintf('%-24s | %-12s | %-12s\n', 'Metric', 'Ege2019', 'Durmaz2024');
fprintf('------------------------------------------------------------\n');
fprintf('%-24s | %-12s | %-12s\n', 'Reached goal', tf2str(results.ege.converged), tf2str(results.durmaz.converged));
fprintf('%-24s | %-12.2f | %-12.2f\n', 'Mission duration (s)', results.ege.duration, results.durmaz.duration);
fprintf('%-24s | %-12.3f | %-12.3f\n', 'Path length (m)', results.ege.path_length, results.durmaz.path_length);
fprintf('%-24s | %-12.3f | %-12.3f\n', 'Average speed (m/s)', results.ege.avg_speed, results.durmaz.avg_speed);
fprintf('%-24s | %-12.4f | %-12.4f\n', 'Avg |yaw rate| (rad/s)', results.ege.avg_abs_yaw_rate, results.durmaz.avg_abs_yaw_rate);
fprintf('%-24s | %-12.4f | %-12.4f\n', 'Final error (m)', results.ege.final_error, results.durmaz.final_error);
fprintf('============================================================\n');

%% PLOT 1: TRAJECTORY OVERLAY ON FUNNEL CHAIN

set(0, 'DefaultFigureVisible', 'on');

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
    plot(ax1, nodes(node_id).poly, 'FaceColor',[1.0 0.85 0.7], 'FaceAlpha',0.20, ...
        'EdgeColor',[1.0 0.5 0.0], 'LineWidth',1.0);
end

h_ege = plot(ax1, results.ege.state_hist(:,1), results.ege.state_hist(:,2), 'r-', 'LineWidth',2.0);
h_durmaz = plot(ax1, results.durmaz.state_hist(:,1), results.durmaz.state_hist(:,2), 'b-', 'LineWidth',2.0);

h_start = plot(ax1, q_start(1), q_start(2), 'go', 'MarkerSize',9, 'LineWidth',2);
text(ax1, q_start(1), q_start(2), '  start', 'FontSize',12, 'FontWeight','bold');

h_goal = plot(ax1, q_goal(1), q_goal(2), 'ko', 'MarkerSize',9, 'LineWidth',2);
text(ax1, q_goal(1), q_goal(2), '  goal', 'FontSize',12, 'FontWeight','bold');

legend([h_ege, h_durmaz, h_start, h_goal], {'Ege2019','Durmaz2024','start','goal'}, 'Location','bestoutside');
title(ax1, 'Trajectory Comparison on Shared Funnel Chain');

%% PLOT 2: SPEED AND YAW RATE

figure('Color','w');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile;
plot(results.ege.time, results.ege.v_hist, 'r-', 'LineWidth',1.5);
hold on;
plot(results.durmaz.time, results.durmaz.v_hist, 'b-', 'LineWidth',1.5);
grid on;
ylabel('v [m/s]');
legend('Ege2019','Durmaz2024','Location','best');
title('Forward Speed');

nexttile;
plot(results.ege.time, results.ege.omega_hist, 'r-', 'LineWidth',1.5);
hold on;
plot(results.durmaz.time, results.durmaz.omega_hist, 'b-', 'LineWidth',1.5);
grid on;
ylabel('\omega [rad/s]');
xlabel('Time [s]');
title('Yaw Rate');

%% FUNCTIONS

function r = packResults(time, state_hist, v_hist, omega_hist, q_goal, goal_tol)

    r.time = time;
    r.state_hist = state_hist;
    r.v_hist = v_hist;
    r.omega_hist = omega_hist;

    r.duration = time(end);
    r.final_error = norm(state_hist(end,1:2) - q_goal);
    r.converged = r.final_error <= goal_tol;

    steps = diff(state_hist(:,1:2), 1, 1);
    r.path_length = sum(hypot(steps(:,1), steps(:,2)));

    r.avg_speed = r.path_length / r.duration;
    r.avg_abs_yaw_rate = trapz(time, abs(omega_hist)) / r.duration;
end

function s = tf2str(tf)
    if tf
        s = 'yes';
    else
        s = 'no';
    end
end

function s = summaryLine(r)
    s = sprintf('%s in %.2fs, final error %.4fm', tf2str(r.converged), r.duration, r.final_error);
end
