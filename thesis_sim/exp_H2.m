%% H2 -- SAFETY FILTERS UNDER UNMODELED DISTURBANCES
%
% H2: Under disturbances that the plan and the controller do not know
%     (current, INS noise, actuator faults), a barrier-function safety
%     filter reduces funnel exits compared with the unfiltered controller.
%
% Controller: radius-scheduled mission speed (MSR, H1), U = 2 m/s (the tug's
% reported top speed), with
%   none  : no filter
%   CBF   : kinematic CBF (surge) / HOCBF (yaw rate) filter on the references
%   HOCBF : second-order HOCBF filter on the thrust (vessel model)
%   rCBF, rHOCBF : the same filters made robust for a current up to
%           cfg.Vb = 0.3 m/s (b_dot replaced by its worst case b_dot - Vb)
% Disturbances (none of them is known to the controller or the filters):
%   nominal   : none
%   current   : 0.15 and 0.3 m/s, 8 directions
%   INS       : low  (0.1 m, 1 deg, 0.02 m/s | rad/s) and
%               high (0.3 m, 3 deg, 0.05) white noise, 5 noise seeds
%   actuator  : left or right thruster at 60 % efficiency, 2 N noise
% Every case is run on every chain; all controllers see exactly the
% same disturbance (same direction / noise seed).
%
% Same map, chains and vessel as H1 and S (thesisConfig.m). Results:
% results/H2_runs.csv, results/H2_summary.csv, figures/H2_*.png.

clear; clc; close all;
cfg = thesisConfig();
rerun = true;                                   % false: replot from results/H2.mat
resFile = fullfile(cfg.dirResults, 'H2.mat');
U = 2.0;
filters = {'none', 'cbf', 'hocbf', 'rcbf', 'rhocbf'};
filterNames = {'none', 'CBF', 'HOCBF', 'rCBF', 'rHOCBF'};
nF = numel(filters);

% Disturbance cases: {level label, dist struct}
cases = {};
cases(end+1, :) = {'nominal', struct()};
for Vc = [0.15 0.3]
    for ang = 0:45:315
        cases(end+1, :) = {sprintf('current %.2f m/s', Vc), struct('Vc', Vc*[cosd(ang); sind(ang)])}; %#ok<SAGROW>
    end
end
ins = {'INS low', 0.1, deg2rad(1), 0.02; 'INS high', 0.3, deg2rad(3), 0.05};
for i = 1:2
    for ns = 1:5
        cases(end+1, :) = {ins{i,1}, struct('ins_pos', ins{i,2}, 'ins_psi', ins{i,3}, 'ins_vel', ins{i,4}, 'seed', ns)}; %#ok<SAGROW>
    end
end
cases(end+1, :) = {'actuator left 60%',  struct('act_gain', [0.6; 1], 'act_noise', 2)};
cases(end+1, :) = {'actuator right 60%', struct('act_gain', [1; 0.6], 'act_noise', 2)};
levels = unique(cases(:,1), 'stable');

if rerun
    jobs = struct('level', {}, 'caseId', {}, 'filter', {}, 'seed', {}, 'run', {});
    for c = 1:size(cases, 1)
        for f = 1:nF
            for s = cfg.seeds
                jobs(end+1) = struct('level', cases{c,1}, 'caseId', c, 'filter', filterNames{f}, 'seed', s, ...
                    'run', struct('law', 'msR', 'U_m', U, 'filter', filters{f}, 'dist', cases{c,2})); %#ok<SAGROW>
            end
        end
    end
    T = runBatch(cfg, jobs);
    save(resFile, 'T', 'cases', 'cfg');
else
    load(resFile, 'T', 'cases');
end

%% ---------------- Summary tables ----------------
T.exited = T.exitAct > 0.01;
T.leftChain = T.exitChain > 0.01;
S = groupsummary(T, {'level', 'filter'}, {'mean', 'max'}, ...
    {'exited', 'exitAct', 'leftChain', 'reached', 'stuck', 'collision', 'T', 'filtPct', 'infeasPct'});
S = S(:, {'level', 'filter', 'GroupCount', 'mean_exited', 'mean_exitAct', 'max_exitAct', 'mean_leftChain', ...
          'mean_reached', 'mean_stuck', 'max_collision', 'mean_T', 'mean_filtPct', 'mean_infeasPct'});
S.Properties.VariableNames = {'level', 'filter', 'runs', 'P_exit', 'mean_exit_m', 'max_exit_m', ...
    'P_left_chain', 'P_reached', 'P_stuck', 'any_collision', 'mean_T_s', 'filter_active_pct', 'infeasible_pct'};
[S.P_exit_lo, S.P_exit_hi] = wilsonCI(S.P_exit, S.runs);          % 95 % intervals
[S.P_reached_lo, S.P_reached_hi] = wilsonCI(S.P_reached, S.runs);
% keep the case order of the script
[~, ord] = sortrows([cellfun(@(l) find(strcmp(levels, l)), S.level), ...
                     cellfun(@(f) find(strcmp(filterNames, f)), S.filter)]);
S = S(ord, :);
writetable(removevars(T, 'funnelExit'), fullfile(cfg.dirResults, 'H2_runs.csv'));
writetable(S, fullfile(cfg.dirResults, 'H2_summary.csv'));
disp(S);

% Paired comparison: same disturbance and chain, filter vs none
P = table();
for f = 2:nF
    for l = 1:numel(levels)
        a = T(strcmp(T.level, levels{l}) & strcmp(T.filter, 'none'), :);
        b = T(strcmp(T.level, levels{l}) & strcmp(T.filter, filterNames{f}), :);
        a = sortrows(a, {'caseId', 'seed'}); b = sortrows(b, {'caseId', 'seed'});
        d = b.exitAct - a.exitAct;
        P = [P; table(levels(l), filterNames(f), sum(d < -0.01), sum(abs(d) <= 0.01), sum(d > 0.01), ...
            'VariableNames', {'level', 'filter', 'better', 'same', 'worse'})]; %#ok<AGROW>
    end
end
writetable(P, fullfile(cfg.dirResults, 'H2_paired.csv'));
disp(P);

%% ---------------- Figures ----------------
cols = lines(nF);
xcat = categorical(levels, levels);

% H2_summary: exit probability, mean exit and goal reached per disturbance
f = figure('Color', 'w', 'Position', [100 100 1500 440]);
tiledlayout(1, 3, 'Padding', 'compact');
st = {'P_exit', 100, 'runs leaving the active funnel [%]'; ...
      'mean_exit_m', 1, 'mean worst exit per run [m]'; ...
      'P_reached', 100, 'runs reaching the goal [%]'};
for p = 1:3
    nexttile; hold on; grid on;
    Y = zeros(numel(levels), nF);
    for l = 1:numel(levels)
        for k = 1:nF
            Y(l, k) = st{p,2}*S.(st{p,1})(strcmp(S.level, levels{l}) & strcmp(S.filter, filterNames{k}));
        end
    end
    b = bar(xcat, Y);
    for k = 1:nF, b(k).FaceColor = cols(k,:); end
    ylabel(st{p,3});
    if p == 1, legend(filterNames, 'Location', 'northwest'); end
end
sgtitle(sprintf('Tugboat, radius-scheduled mission speed U = %.0f m/s, %d chains', U, numel(cfg.seeds)));
saveFigure(f, cfg, 'H2_summary');

% H2_activity: how often the filters act and how often they cannot satisfy the barrier
f = figure('Color', 'w', 'Position', [100 100 1500 420]);
tiledlayout(1, 3, 'Padding', 'compact');
st = {'filter_active_pct', 'steps the filter changed the command [%]'; ...
      'infeasible_pct',    'steps the barrier condition was unsatisfiable [%]'; ...
      'P_stuck',           'share of runs stalled (no progress for 300 s)'};
for p = 1:3
    nexttile; hold on; grid on;
    Y = zeros(numel(levels), nF-1);
    for l = 1:numel(levels)
        for k = 2:nF
            Y(l, k-1) = S.(st{p,1})(strcmp(S.level, levels{l}) & strcmp(S.filter, filterNames{k}));
        end
    end
    b = bar(xcat, Y);
    for k = 1:nF-1, b(k).FaceColor = cols(k+1,:); end
    ylabel(st{p,2});
    if p == 1, legend(filterNames(2:nF), 'Location', 'northwest'); end
end
saveFigure(f, cfg, 'H2_activity');

% H2_example: the disturbance case with the largest unfiltered exit
Tn = T(strcmp(T.filter, 'none'), :);
[~, iw] = max(Tn.exitAct);
cw = Tn.caseId(iw); sw = Tn.seed(iw); ch = buildChain(cfg, sw);
f = figure('Color', 'w', 'Position', [100 100 1000 380]); hold on; grid on;
for k = 1:nF
    o = simulateRun(ch, cfg, struct('law', 'msR', 'U_m', U, 'filter', filters{k}, ...
        'dist', cases{cw, 2}, 'keepLog', true));
    plot(o.log.t, o.log.h, 'Color', cols(k,:), 'LineWidth', 1.2);
end
yline(0, 'k--');
xlabel('t [s]'); ylabel('h = R_{active} - \rho [m]  (h < 0: outside)');
legend(filterNames, 'Location', 'best');
title(sprintf('Worst unfiltered case: %s, chain %d', cases{cw, 1}, sw));
saveFigure(f, cfg, 'H2_example');
