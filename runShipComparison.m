function runShipComparison()
%RUNSHIPCOMPARISON  Durmaz vs CBF vs HOCBF on the RSC map for three ships.
%
% Runs compareFunnelChain.m once per ship with
%   mission speed U_m = the ship's maximum speed (vessel.Umax),
%   radius speed scaling on (U_k = min(U_m, R_k*Ka)), tanh goal profile,
% and saves its figures to results/<ship>_<figure>.png plus a summary
% table to results/shipComparison.mat. The map is whatever RSC.m builds
% (scenarioId in RSC.m); all ships use the same funnel chain.

    ships  = {'tugboat', 'otter', 'cybership'};
    outDir = fullfile(fileparts(mfilename('fullpath')), 'results');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    rows = {};
    for i = 1:numel(ships)
        p = feval([ships{i} '3d']);
        cfg = struct('vessel', ships{i}, 'U_m', p.Umax, 'profile', 'tanh', 'scaleR', true);
        fprintf('\n===== %s: U_m = %.2f m/s =====\n', ships{i}, p.Umax);

        % compareFunnelChain.m is a script that clears the base workspace,
        % so it runs there with cfg as its only input
        assignin('base', 'cfg', cfg);
        evalin('base', 'compareFunnelChain');

        % Save every figure it made
        figs = findall(0, 'Type', 'figure');
        for f = figs'
            name = lower(regexprep(f.Name, '^Funnel chain - ', ''));
            name = regexprep(name, '[^a-z0-9]+', '_');
            exportgraphics(f, fullfile(outDir, sprintf('%s_%s.png', ships{i}, name)), 'Resolution', 200);
        end
        close(figs);

        % Collect the results table
        res = evalin('base', 'res'); ctrls = evalin('base', 'ctrls');
        for c = 1:numel(ctrls)
            r = res{c};
            rows(end+1, :) = {ships{i}, p.Umax, ctrls{c}, r.reached, r.time(end), ...
                r.viol, r.t_out, r.path_len}; %#ok<AGROW>
        end
    end

    summary = cell2table(rows, 'VariableNames', ...
        {'ship', 'U_m', 'controller', 'reached', 'T', 'max_viol', 't_out', 'path_len'});
    save(fullfile(outDir, 'shipComparison.mat'), 'summary');
    fprintf('\n===== Summary =====\n');
    disp(summary);
    fprintf('Plots and summary saved to %s\n', outDir);
end
