%% H2 -- RADIUS-SCHEDULED MISSION SPEED
%
% H2: With a constant mission speed U the tugboat leaves funnels whose radius
%     is small compared with the distance it travels while turning
%     (R_k/U below the heading time constant T_R = 1/Ka). Scheduling the
%     cruise speed by funnel radius, U_k = min(U, R_k/T_R), removes these
%     exits at a small cost in travel time.
%
% Runs: 'ms' (constant U) and 'msR' (radius-scheduled) at U = 1 ... 3.5 m/s
% on every chain, no disturbance. U > 2 m/s is above the tug's reported top
% speed (Erunsal et al. 2017) but within the model's 3.7 m/s at full thrust;
% it is included to show where the exits start.
%
% Same map, chains and vessel as H1 and H3 (thesisConfig.m). Results:
% results/H2_runs.csv, results/H2_summary.csv, figures/H2_*.png.

clear; clc; close all;
cfg = thesisConfig();
rerun = true;                                   % false: replot from results/H2.mat
resFile = fullfile(cfg.dirResults, 'H2.mat');
Us = [1 1.5 2 2.5 3 3.5];
lawIds = {'ms', 'msR'}; lawNames = {'MS (constant U)', 'MS + radius scheduling'};

if rerun
    jobs = struct('law', {}, 'U', {}, 'seed', {}, 'run', {});
    for i = 1:2
        for U = Us
            for s = cfg.seeds
                jobs(end+1) = struct('law', lawNames{i}, 'U', U, 'seed', s, ...
                    'run', struct('law', lawIds{i}, 'U_m', U)); %#ok<SAGROW>
            end
        end
    end
    T = runBatch(cfg, jobs);
    chains = arrayfun(@(s) buildChain(cfg, s), cfg.seeds, 'UniformOutput', false);
    save(resFile, 'T', 'chains', 'cfg');
else
    load(resFile, 'T', 'chains');
end

%% ---------------- Summary tables ----------------
T.exited = T.exitAct > 0.01;
S = groupsummary(T, {'law', 'U'}, {'mean', 'max'}, {'exited', 'exitAct', 'exitChain', 'reached', 'collision', 'T'});
S = S(:, {'law', 'U', 'GroupCount', 'mean_exited', 'max_exitAct', 'max_exitChain', 'mean_reached', ...
          'max_collision', 'mean_T'});
S.Properties.VariableNames = {'law', 'U', 'runs', 'P_exit', 'max_exit_m', 'max_exit_chain_m', ...
                              'P_reached', 'any_collision', 'mean_T_s'};
writetable(removevars(T, 'funnelExit'), fullfile(cfg.dirResults, 'H2_runs.csv'));
writetable(S, fullfile(cfg.dirResults, 'H2_summary.csv'));
disp(S);

%% ---------------- Figures ----------------
cols = lines(2);

% H2_sweep: exit probability, worst exit and travel time vs mission speed
f = figure('Color', 'w', 'Position', [100 100 1200 360]);
tiledlayout(1, 3, 'Padding', 'compact');
st = {@(r) 100*mean(r.exited), 'runs leaving a funnel [%]'; ...
      @(r) max(r.exitAct),     'worst exit from the active funnel [m]'; ...
      @(r) mean(r.T),          'mean time to goal [s]'};
for p = 1:3
    nexttile; hold on; grid on;
    for i = 1:2
        y = arrayfun(@(U) st{p,1}(T(strcmp(T.law, lawNames{i}) & T.U == U, :)), Us);
        plot(Us, y, '-o', 'Color', cols(i,:), 'MarkerFaceColor', cols(i,:), 'LineWidth', 1.5);
    end
    xline(2, 'k:', 'reported top speed', 'LabelVerticalAlignment', 'bottom');
    xlabel('mission speed U [m/s]'); ylabel(st{p,2});
    if p == 1, legend(lawNames, 'Location', 'northwest'); end
end
sgtitle(sprintf('Tugboat, no disturbance (%d chains per point)', numel(cfg.seeds)));
saveFigure(f, cfg, 'H2_sweep');

% H2_mechanism: exit from each funnel vs R_k / (speed in that funnel)
f = figure('Color', 'w', 'Position', [100 100 900 380]);
tiledlayout(1, 2, 'Padding', 'compact');
for i = 1:2
    nexttile; hold on; grid on;
    rows = T(strcmp(T.law, lawNames{i}), :);
    x = []; y = [];
    for j = 1:height(rows)
        R = chains{cfg.seeds == rows.seed(j)}.R;
        Uk = rows.U(j)*ones(size(R));
        if i == 2, Uk = min(rows.U(j), R/cfg.T_R); end
        e = rows.funnelExit{j};
        ok = ~isnan(e); ok(end) = false;          % visited intermediate funnels
        x = [x, R(ok)./Uk(ok)]; y = [y, e(ok)]; %#ok<AGROW>
    end
    scatter(x, y, 14, cols(i,:), 'filled', 'MarkerFaceAlpha', 0.5);
    xline(cfg.T_R, 'k--', 'T_R = 1/K_a', 'LabelVerticalAlignment', 'top');
    set(gca, 'XScale', 'log');
    xlabel('R_k / U_k  [s]'); ylabel('exit from funnel k [m]');
    title(lawNames{i});
end
sgtitle('Each point: one funnel in one run (all speeds, all chains)');
saveFigure(f, cfg, 'H2_mechanism');

% H2_example: trajectories at the highest speed on the chain with the smallest funnel
[~, iMin] = min(cellfun(@(c) min(c.R), chains));
ch = chains{iMin}; Umax = Us(end);
f = figure('Color', 'w', 'Position', [100 100 700 650]); hold on; axis equal; box on;
th = linspace(0, 2*pi, 100);
for k = 1:numel(ch.obs), plot(ch.obs{k}, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none'); end
for k = 1:numel(ch.R)
    plot(ch.C(1,k) + ch.R(k)*cos(th), ch.C(2,k) + ch.R(k)*sin(th), 'Color', [1 0.6 0.2]);
end
h = gobjects(1, 2);
for i = 1:2
    o = simulateRun(ch, cfg, struct('law', lawIds{i}, 'U_m', Umax, 'keepLog', true));
    h(i) = plot(o.log.X, o.log.Y, 'Color', cols(i,:), 'LineWidth', 1.3);
end
plot(ch.q_start(1), ch.q_start(2), 'go', 'MarkerFaceColor', 'g');
plot(ch.q_goal(1), ch.q_goal(2), 'ro', 'MarkerFaceColor', 'r');
xlim(ch.W(1:2)); ylim(ch.W(3:4));
legend(h, lawNames, 'Location', 'southoutside');
title(sprintf('Chain %d (smallest funnel %.1f m), U = %.1f m/s', ch.seed, min(ch.R), Umax));
saveFigure(f, cfg, 'H2_example');
