%% S -- SENSITIVITY OF THE RADIUS-SCHEDULED MISSION SPEED (MSR)
%
% Do MSR's results (no funnel exits, travel time) depend on
%   1. the scheduling constant T_R (U_k = min(U, R_k/T_R)), here 0.5, 1 and 2
%      times the heading time constant 1/Ka;
%   2. the tugboat's reverse thrust, which was not measured (assumed -14.5 N),
%      here -7.25, -14.5 and -26 N;
%   3. the thruster time constant, which was not identified (assumed 0.25 s),
%      here 0.1, 0.25 and 0.5 s?
% Each variant: MSR at U = 2 and 3 m/s (the tug's reported top speed and
% above it), no disturbance, every chain.
%
% Results: results/S_runs.csv, results/S_summary.csv, figures/S_sensitivity.png.

clear; clc; close all;
cfg = thesisConfig();
rerun = true;                                   % false: replot from results/S.mat
resFile = fullfile(cfg.dirResults, 'S.mat');
Us = [2 3];

% {part, label, value, config, laws}
V = {};
for m = [0.5 1 2]
    c = thesisConfig(); c.T_R = m/c.Ka;
    V(end+1, :) = {'T_R', sprintf('T_R = %.1f/K_a', m), m, c, {'msR'}}; %#ok<SAGROW>
end
for Fmin = [-7.25 -14.5 -26]
    V(end+1, :) = {'Fmin', sprintf('F_{min} = %.1f N', Fmin), Fmin, ...
        thesisConfig('tugboat', struct('Fmin', Fmin)), {'msR'}}; %#ok<SAGROW>
end
for Ta = [0.1 0.25 0.5]
    V(end+1, :) = {'T_act', sprintf('T_{act} = %.2f s', Ta), Ta, ...
        thesisConfig('tugboat', struct('T_act', Ta)), {'msR'}}; %#ok<SAGROW>
end
lawName = containers.Map({'ms', 'msR'}, {'MS (constant U)', 'MS + radius scheduling'});

if rerun
    T = table();
    for v = 1:size(V, 1)
        jobs = struct('part', {}, 'value', {}, 'law', {}, 'U', {}, 'seed', {}, 'run', {});
        for l = V{v,5}
            for U = Us
                for s = cfg.seeds
                    jobs(end+1) = struct('part', V{v,1}, 'value', V{v,3}, 'law', lawName(l{1}), ...
                        'U', U, 'seed', s, 'run', struct('law', l{1}, 'U_m', U)); %#ok<SAGROW>
                end
            end
        end
        fprintf('%s\n', V{v,2});
        T = [T; runBatch(V{v,4}, jobs)]; %#ok<AGROW>
    end
    save(resFile, 'T');
else
    load(resFile, 'T');
end

%% ---------------- Summary table ----------------
T.exited = T.exitAct > 0.01;
S = groupsummary(T, {'part', 'value', 'law', 'U'}, {'mean', 'max'}, {'exited', 'exitAct', 'reached', 'T'});
S = S(:, {'part', 'value', 'law', 'U', 'GroupCount', 'mean_exited', 'max_exitAct', 'mean_reached', 'mean_T'});
S.Properties.VariableNames = {'part', 'value', 'law', 'U', 'runs', 'P_exit', 'max_exit_m', 'P_reached', 'mean_T_s'};
[S.P_exit_lo, S.P_exit_hi] = wilsonCI(S.P_exit, S.runs);
writetable(removevars(T, 'funnelExit'), fullfile(cfg.dirResults, 'S_runs.csv'));
writetable(S, fullfile(cfg.dirResults, 'S_summary.csv'));
disp(S);

%% ---------------- Figure ----------------
% Top: share of runs leaving a funnel; bottom: mean time to goal
parts = {'T_R', 'Fmin', 'T_act'};
xl = {'T_R \cdot K_a  (scheduling constant / heading time constant)', ...
      'assumed reverse thrust F_{min} [N]', 'assumed thruster time constant [s]'};
cols = lines(2); ln = lawName('msR');
f = figure('Color', 'w', 'Position', [100 100 1400 650]);
tiledlayout(2, 3, 'Padding', 'compact');
for row = 1:2
    for p = 1:3
        nexttile; hold on; grid on;
        Sp = S(strcmp(S.part, parts{p}) & strcmp(S.law, ln), :);
        for j = 1:numel(Us)
            r = sortrows(Sp(Sp.U == Us(j), :), 'value');
            if row == 1
                errorbar(r.value, 100*r.P_exit, 100*(r.P_exit - r.P_exit_lo), 100*(r.P_exit_hi - r.P_exit), ...
                    '-o', 'Color', cols(j,:), 'MarkerFaceColor', cols(j,:), 'LineWidth', 1.3);
            else
                plot(r.value, r.mean_T_s, '-o', 'Color', cols(j,:), 'MarkerFaceColor', cols(j,:), 'LineWidth', 1.3);
            end
        end
        xlabel(xl{p});
        if row == 1
            ylabel('runs leaving a funnel [%] (95 % CI)'); ylim([-5 100]);
        else
            ylabel('mean time to goal [s]');
        end
        if row == 1 && p == 1
            legend(arrayfun(@(U) sprintf('MSR, U = %.0f m/s', U), Us, 'UniformOutput', false), 'Location', 'best');
        end
    end
end
sgtitle(sprintf('Tugboat, radius-scheduled mission speed, no disturbance, %d chains per point', numel(cfg.seeds)));
saveFigure(f, cfg, 'S_sensitivity');
