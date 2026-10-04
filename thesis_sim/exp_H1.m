%% H1 -- RADIUS-SCHEDULED MISSION SPEED vs THE ORIGINAL DURMAZ2024 SPEED LAW
%
% H1: Replacing Durmaz2024's speed s = 2*Kv*rho by a mission speed scheduled
%     by funnel radius, U_k = min(U, R_k/T_R), keeps the vessel inside the
%     funnels like the original law does, but removes its stop-and-go speed
%     profile, shortens the travel time and lets the vessel reach the goal
%     against a current.
%
% A. Kinematic guarantee: ideal unicycle (references realized exactly),
%    Durmaz vs MSR at 1, 2, 3 m/s -> exits must be 0 for every law.
% B. Tugboat, no disturbance: Durmaz vs MSR at U = 1, 1.5, 2 m/s (the tug's
%    reported top speed): exits, travel time, speed oscillation
%    CV(u) = std(u)/mean(u) in the intermediate funnels.
% C. Tugboat against a current (0.15 and 0.3 m/s, 8 directions): Durmaz vs
%    MSR at 2 m/s: goal reached, exits.
% D. Otter, no disturbance: Durmaz vs MSR at U = 1, 2, 3 m/s.
%
% Same map, chains and controller as H2 and S (thesisConfig.m). Results:
% results/H1_runs.csv, results/H1_summary.csv, figures/H1_*.png.

clear; clc; close all;
cfg = thesisConfig();
rerun = true;                                   % false: replot from results/H1.mat
resFile = fullfile(cfg.dirResults, 'H1.mat');

durmaz = struct('name', 'Durmaz', 'run', struct('law', 'durmaz'));
msr = @(U) struct('name', sprintf('MSR %.1f m/s', U), 'run', struct('law', 'msR', 'U_m', U));
lawsUni = [durmaz, msr(1), msr(2), msr(3)];
lawsTug = [durmaz, msr(1), msr(1.5), msr(2)];
lawsOtt = [durmaz, msr(1), msr(2), msr(3)];

if rerun
    jobs = struct('part', {}, 'vessel', {}, 'law', {}, 'Vc', {}, 'ang', {}, 'seed', {}, 'run', {});
    add = @(jobs, part, vessel, l, Vc, ang, s, r) [jobs, struct('part', part, 'vessel', vessel, ...
        'law', l.name, 'Vc', Vc, 'ang', ang, 'seed', s, 'run', r)];
    for s = cfg.seeds
        for l = lawsUni                         % A. unicycle
            r = l.run; r.plant = 'unicycle';
            jobs = add(jobs, 'A', 'unicycle', l, 0, 0, s, r);
        end
        for l = lawsTug                         % B. tugboat, nominal
            jobs = add(jobs, 'B', 'tugboat', l, 0, 0, s, l.run);
        end
        for l = [durmaz, msr(2)]                % C. tugboat, current
            for Vc = [0.15 0.3]
                for ang = 0:45:315
                    r = l.run; r.dist = struct('Vc', Vc*[cosd(ang); sind(ang)]);
                    jobs = add(jobs, 'C', 'tugboat', l, Vc, ang, s, r);
                end
            end
        end
    end
    T = runBatch(cfg, jobs);

    % D. Otter, nominal (own vessel configuration)
    cfgO = thesisConfig('otter');
    jobsO = struct('part', {}, 'vessel', {}, 'law', {}, 'Vc', {}, 'ang', {}, 'seed', {}, 'run', {});
    for s = cfg.seeds
        for l = lawsOtt
            jobsO = add(jobsO, 'D', 'otter', l, 0, 0, s, l.run);
        end
    end
    T = [T; runBatch(cfgO, jobsO)];

    % Example time histories (tugboat, first chain, no disturbance)
    ch = buildChain(cfg, cfg.seeds(1));
    EX = struct('name', {}, 'log', {});
    for l = [durmaz, msr(2)]
        r = l.run; r.keepLog = true;
        o = simulateRun(ch, cfg, r);
        EX(end+1) = struct('name', l.name, 'log', o.log); %#ok<SAGROW>
    end
    save(resFile, 'T', 'EX', 'ch');
else
    load(resFile, 'T', 'EX', 'ch');
end

%% ---------------- Summary tables ----------------
T.exited = T.exitAct > 0.01;
S = groupsummary(T, {'part', 'vessel', 'law', 'Vc'}, {'mean', 'max'}, ...
    {'exited', 'exitAct', 'reached', 'collision', 'T', 'uMean', 'uCV'});
S = S(:, {'part', 'vessel', 'law', 'Vc', 'GroupCount', 'mean_exited', 'max_exitAct', 'mean_reached', ...
          'max_collision', 'mean_T', 'mean_uMean', 'mean_uCV'});
S.Properties.VariableNames = {'part', 'vessel', 'law', 'Vc', 'runs', 'P_exit', 'max_exit_m', 'P_reached', ...
                              'any_collision', 'mean_T_s', 'mean_u', 'mean_CV_u'};
[S.P_exit_lo, S.P_exit_hi] = wilsonCI(S.P_exit, S.runs);          % 95 % intervals
[S.P_reached_lo, S.P_reached_hi] = wilsonCI(S.P_reached, S.runs);
writetable(removevars(T, 'funnelExit'), fullfile(cfg.dirResults, 'H1_runs.csv'));
writetable(S, fullfile(cfg.dirResults, 'H1_summary.csv'));
disp(S);

%% ---------------- Figures ----------------
cD = [0 0.45 0.74]; cM = [0.85 0.33 0.10];      % Durmaz, MSR

% H1_trajectory: both laws on the harbour map (tugboat, first chain)
f = figure('Color', 'w', 'Position', [100 100 800 560]); hold on; axis equal; box on;
th = linspace(0, 2*pi, 100);
for k = 1:numel(ch.obs), plot(ch.obs{k}, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none'); end
for k = 1:numel(ch.R)
    plot(ch.C(1,k) + ch.R(k)*cos(th), ch.C(2,k) + ch.R(k)*sin(th), 'Color', [1 0.6 0.2]);
end
h = gobjects(1, 2); cc = [cD; cM];
for i = 1:2
    h(i) = plot(EX(i).log.X, EX(i).log.Y, 'Color', cc(i,:), 'LineWidth', 1.4);
end
plot(ch.q_start(1), ch.q_start(2), 'go', 'MarkerFaceColor', 'g');
plot(ch.q_goal(1), ch.q_goal(2), 'ro', 'MarkerFaceColor', 'r');
xlim(ch.W(1:2)); ylim(ch.W(3:4)); xlabel('x [m]'); ylabel('y [m]');
legend(h, {EX.name}, 'Location', 'southeast');
title(sprintf('Tugboat, chain %d', ch.seed));
saveFigure(f, cfg, 'H1_trajectory');

% H1_speed: surge speed over time
f = figure('Color', 'w', 'Position', [100 100 900 340]); hold on; grid on;
for i = 1:2
    plot(EX(i).log.t, EX(i).log.u, 'Color', cc(i,:), 'LineWidth', 1.2);
end
xlabel('t [s]'); ylabel('surge speed u [m/s]');
legend({EX.name}, 'Location', 'northeast');
title(sprintf('Tugboat, chain %d: original speed law vs radius-scheduled mission speed', ch.seed));
saveFigure(f, cfg, 'H1_speed');

% H1_nominal: travel time, speed oscillation and worst exit per law (tugboat and Otter)
f = figure('Color', 'w', 'Position', [100 100 1300 680]);
tiledlayout(2, 3, 'Padding', 'compact');
parts = {'B', 'D'}; laws = {lawsTug, lawsOtt}; vn = {'Tugboat', 'Otter'};
st = {'T', 'time to goal [s]'; 'uCV', 'CV(u) = std(u)/mean(u)'; 'exitAct', 'worst exit [m]'};
for v = 1:2
    names = {laws{v}.name};
    Tv = T(strcmp(T.part, parts{v}), :);
    for p = 1:3
        nexttile; hold on; grid on;
        y = cellfun(@(n) mean(Tv.(st{p,1})(strcmp(Tv.law, n))), names);
        if p == 3, y = cellfun(@(n) max(Tv.(st{p,1})(strcmp(Tv.law, n))), names); end
        short = regexprep(names, '\.0? m/s| m/s', '');   % 'MSR 1.5 m/s' -> 'MSR 1.5'
        b = bar(categorical(short, short), y, 'FaceColor', 'flat');
        xtickangle(0); xlabel('MSR: mission speed U [m/s]');
        b.CData = [cD; repmat(cM, numel(names) - 1, 1)];
        ylabel(st{p,2}); title(vn{v});
        if p == 3
            ylim([0 max(0.1, 1.2*max(y))]);
            if all(y == 0), text(2.5, 0.05, 'no exits', 'HorizontalAlignment', 'center', 'FontSize', 12); end
        end
    end
end
sgtitle(sprintf('No disturbance, mean over %d chains (worst exit: max)', numel(cfg.seeds)));
saveFigure(f, cfg, 'H1_nominal');

% H1_current: goal reached and exits against a current (tugboat)
f = figure('Color', 'w', 'Position', [100 100 900 360]);
tiledlayout(1, 2, 'Padding', 'compact');
Cc = T(strcmp(T.part, 'C'), :); Vs = unique(Cc.Vc)';
m2 = msr(2); names = {durmaz.name, m2.name}; cc = [cD; cM];
for p = 1:2
    nexttile; hold on; grid on;
    for i = 1:2
        sel = @(v) strcmp(Cc.law, names{i}) & Cc.Vc == v;
        if p == 1
            pr = arrayfun(@(v) mean(Cc.reached(sel(v))), Vs); nr = arrayfun(@(v) sum(sel(v)), Vs);
            [lo, hi] = wilsonCI(pr, nr);
            errorbar(Vs, 100*pr, 100*(pr - lo), 100*(hi - pr), '-o', 'Color', cc(i,:), ...
                'MarkerFaceColor', cc(i,:), 'LineWidth', 1.5);
        else
            plot(Vs, arrayfun(@(v) max(Cc.exitAct(sel(v))), Vs), '-o', 'Color', cc(i,:), ...
                'MarkerFaceColor', cc(i,:), 'LineWidth', 1.5);
        end
    end
    xlim([min(Vs) - 0.05, max(Vs) + 0.05]); xlabel('current speed [m/s]');
    if p == 1
        ylabel('runs reaching the goal [%] (95 % CI)'); ylim([-5 105]);
        legend(names, 'Location', 'southwest');
    else
        ylabel('worst exit [m]');
    end
end
sgtitle(sprintf('Tugboat against a current (8 directions x %d chains)', numel(cfg.seeds)));
saveFigure(f, cfg, 'H1_current');
