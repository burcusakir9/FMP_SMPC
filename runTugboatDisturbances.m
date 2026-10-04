function runTugboatDisturbances()
%RUNTUGBOATDISTURBANCES  Durmaz vs CBF vs HOCBF on the tugboat under disturbances.
%
% Runs compareFunnelChain.m for the tugboat (U_m = 2 m/s, radius speed
% scaling, tanh goal profile) once per disturbance case below, and saves
% its figures to results/tugboat_<case>_<figure>.png plus a summary table
% to results/tugboatDisturbances.mat. None of the disturbances is known to
% the controller or the filters (see the Disturbances section of
% compareFunnelChain.m).

    cases = struct( ...
        'name', {'nominal', 'current', 'ins', 'actuator'}, ...
        'dist', {struct(), ...
                 struct('Vc', [0; -0.3]), ...                                  % 0.3 m/s current
                 struct('ins_pos', 0.3, 'ins_psi', deg2rad(3), 'ins_vel', 0.05), ...
                 struct('act_gain', [0.6; 1.0], 'act_noise', 2)});              % weak left motor + noise

    outDir = fullfile(fileparts(mfilename('fullpath')), 'results');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    rows = {};
    for i = 1:numel(cases)
        cfg = struct('vessel', 'tugboat', 'U_m', 2.0, 'profile', 'tanh', 'scaleR', true, ...
            'dist', cases(i).dist);
        fprintf('\n===== tugboat, disturbance: %s =====\n', cases(i).name);

        % compareFunnelChain.m is a script that clears the base workspace,
        % so it runs there with cfg as its only input
        assignin('base', 'cfg', cfg);
        evalin('base', 'compareFunnelChain');

        % Save every figure it made
        figs = findall(0, 'Type', 'figure');
        for f = figs'
            name = lower(regexprep(f.Name, '^Funnel chain - ', ''));
            name = regexprep(name, '[^a-z0-9]+', '_');
            exportgraphics(f, fullfile(outDir, sprintf('tugboat_%s_%s.png', cases(i).name, name)), ...
                'Resolution', 200);
        end
        close(figs);

        % Collect the results table
        res = evalin('base', 'res'); ctrls = evalin('base', 'ctrls');
        for c = 1:numel(ctrls)
            r = res{c};
            rows(end+1, :) = {cases(i).name, ctrls{c}, r.reached, r.time(end), ...
                r.viol, r.t_out, 100*r.filt_frac}; %#ok<AGROW>
        end
    end

    summary = cell2table(rows, 'VariableNames', ...
        {'disturbance', 'controller', 'reached', 'T', 'max_viol', 't_out', 'filter_pct'});
    save(fullfile(outDir, 'tugboatDisturbances.mat'), 'summary', 'cases');
    fprintf('\n===== Summary =====\n');
    disp(summary);
    fprintf('Plots and summary saved to %s\n', outDir);
end
