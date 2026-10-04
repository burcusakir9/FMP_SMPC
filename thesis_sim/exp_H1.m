%% H1 -- MISSION SPEED vs THE ORIGINAL DURMAZ2024 SPEED LAW
%
% H1: Replacing Durmaz2024's speed s = 2*Kv*rho by a mission (cruise) speed
%     keeps the kinematic funnel guarantee, and on the tugboat removes the
%     stop-and-go speed profile and shortens the travel time.
%
% A. Kinematic guarantee: ideal unicycle (references realized exactly),
%    Durmaz vs mission speed at 1, 2, 3 m/s -> exits must be 0 for every law.
% B. Tugboat dynamics, no disturbance: travel time, mean speed and speed
%    oscillation CV(u) = std(u)/mean(u) in the intermediate funnels.
% C. Robustness to a current (0.15 and 0.3 m/s, 8 directions): goal reached?
%    The original law slows down near every center, where the current wins.
%
% Same map, chains and vessel as H2 and H3 (thesisConfig.m). Results:
% results/H1_runs.csv, results/H1_summary.csv, figures/H1_*.png.

clear; clc; close all;
cfg = thesisConfig();
rerun = true;                                   % false: replot from results/H1.mat
resFile = fullfile(cfg.dirResults, 'H1.mat');

laws = struct('name', {'Durmaz', 'MS 1 m/s', 'MS 2 m/s'}, ...
              'run',  {struct('law', 'durmaz'), struct('law', 'ms', 'U_m', 1), struct('law', 'ms', 'U_m', 2)});

if rerun
    jobs = struct('part', {}, 'law', {}, 'Vc', {}, 'ang', {}, 'seed', {}, 'run', {});
    % A. unicycle
    uniLaws = [laws, struct('name', 'MS 3 m/s', 'run', struct('law', 'ms', 'U_m', 3))];
    for l = uniLaws
        for s = cfg.seeds
            r = l.run; r.plant = 'unicycle';
            jobs(end+1) = struct('part', 'A', 'law', l.name, 'Vc', 0, 'ang', 0, 'seed', s, 'run', r); %#ok<SAGROW>
        end
    end
    % B. tugboat, nominal
    for l = laws
        for s = cfg.seeds
            jobs(end+1) = struct('part', 'B', 'law', l.name, 'Vc', 0, 'ang', 0, 'seed', s, 'run', l.run); %#ok<SAGROW>
        end
    end
    % C. tugboat, current
    for l = laws([1 3])
        for Vc = [0.15 0.3]
            for ang = 0:45:315
                for s = cfg.seeds
                    r = l.run; r.dist = struct('Vc', Vc*[cosd(ang); sind(ang)]);
                    jobs(end+1) = struct('part', 'C', 'law', l.name, 'Vc', Vc, 'ang', ang, 'seed', s, 'run', r); %#ok<SAGROW>
                end
            end
        end
    end
    T = runBatch(cfg, jobs);

    % Example time histories (chain of the first seed, no disturbance)
    ch = buildChain(cfg, cfg.seeds(1));
    EX = struct('name', {}, 'log', {});
    for l = laws([1 3])
        r = l.run; r.keepLog = true;
        o = simulateRun(ch, cfg, r);
        EX(end+1) = struct('name', l.name, 'log', o.log); %#ok<SAGROW>
    end
    save(resFile, 'T', 'EX', 'cfg');
else
    load(resFile, 'T', 'EX');
end

%% ---------------- Summary tables ----------------
T.exited = T.exitAct > 0.01;
S = groupsummary(T, {'part', 'law', 'Vc'}, {'mean', 'max'}, ...
    {'exited', 'exitAct', 'reached', 'collision', 'T', 'uMean', 'uCV'});
S = S(:, {'part', 'law', 'Vc', 'GroupCount', 'mean_exited', 'max_exitAct', 'mean_reached', ...
          'max_collision', 'mean_T', 'mean_uMean', 'mean_uCV'});
S.Properties.VariableNames = {'part', 'law', 'Vc', 'runs', 'P_exit', 'max_exit_m', 'P_reached', ...
                              'any_collision', 'mean_T_s', 'mean_u', 'mean_CV_u'};
writetable(removevars(T, 'funnelExit'), fullfile(cfg.dirResults, 'H1_runs.csv'));
writetable(S, fullfile(cfg.dirResults, 'H1_summary.csv'));
disp(S);

%% ---------------- Figures ----------------
cols = lines(3); lawNames = {laws.name};

% H1_speed: example surge speed, original law vs mission speed
f = figure('Color', 'w', 'Position', [100 100 900 360]); hold on; grid on;
for i = 1:numel(EX)
    plot(EX(i).log.t, EX(i).log.u, 'Color', cols(2*i-1,:), 'LineWidth', 1.2);
end
xlabel('t [s]'); ylabel('surge speed u [m/s]');
legend({EX.name}, 'Location', 'best');
title(sprintf('Tugboat, chain %d: original speed law vs mission speed', cfg.seeds(1)));
saveFigure(f, cfg, 'H1_speed');

% H1_nominal: travel time and speed oscillation per law (mean +- std over chains)
f = figure('Color', 'w', 'Position', [100 100 900 360]);
tiledlayout(1, 2, 'Padding', 'compact');
B = T(strcmp(T.part, 'B'), :);
for p = 1:2
    nexttile; hold on; grid on;
    v = {'T', 'uCV'}; lab = {'time to goal [s]', 'CV(u) = std(u)/mean(u)'};
    m = cellfun(@(n) mean(B.(v{p})(strcmp(B.law, n))), lawNames);
    s = cellfun(@(n) std(B.(v{p})(strcmp(B.law, n))), lawNames);
    bar(categorical(lawNames, lawNames), m, 'FaceColor', [0.3 0.5 0.8]);
    errorbar(1:3, m, s, 'k.', 'LineWidth', 1);
    ylabel(lab{p});
end
sgtitle(sprintf('Tugboat, no disturbance (%d chains)', numel(cfg.seeds)));
saveFigure(f, cfg, 'H1_nominal');

% H1_current: share of runs reaching the goal against a current
f = figure('Color', 'w', 'Position', [100 100 500 360]); hold on; grid on;
Cc = T(strcmp(T.part, 'C'), :); Vs = unique(Cc.Vc)';
for i = [1 3]
    y = arrayfun(@(v) 100*mean(Cc.reached(strcmp(Cc.law, lawNames{i}) & Cc.Vc == v)), Vs);
    plot(Vs, y, '-o', 'Color', cols(i,:), 'MarkerFaceColor', cols(i,:), 'LineWidth', 1.5);
end
xlabel('current speed [m/s]'); ylabel('runs reaching the goal [%]'); ylim([-5 105]);
legend(lawNames([1 3]), 'Location', 'southwest');
title('Tugboat against a current (8 directions x chains)');
saveFigure(f, cfg, 'H1_current');
