function T = runBatch(cfg, jobs)
%RUNBATCH  Run a list of simulations (in parallel if a pool is available).
%
%   T = runBatch(cfg, jobs)
%
% jobs: struct array, each with fields
%   seed   RSC chain seed (buildChain)
%   run    run struct for simulateRun
%   plus any label fields (e.g. case, level, filter), copied to the table.
% T: one table row per job with the labels and the simulateRun metrics.

    chains = arrayfun(@(s) buildChain(cfg, s), cfg.seeds, 'UniformOutput', false);
    seedIdx = arrayfun(@(j) find(cfg.seeds == j.seed), jobs);
    n = numel(jobs);
    M = cell(n, 1);
    fprintf('Running %d simulations ...\n', n);
    parfor i = 1:n
        o = simulateRun(chains{seedIdx(i)}, cfg, jobs(i).run);
        M{i} = rmfield(o, 'funnelExit');
        M{i}.funnelExit = {o.funnelExit};          % per-funnel exits (cell, for the table)
    end
    T = struct2table([M{:}]');

    % Label columns first
    labels = setdiff(fieldnames(jobs), {'run'}, 'stable');
    for k = numel(labels):-1:1
        T = addvars(T, {jobs.(labels{k})}', 'Before', 1, 'NewVariableNames', labels{k});
        if all(cellfun(@(v) isnumeric(v) && isscalar(v), T.(labels{k})))
            T.(labels{k}) = cell2mat(T.(labels{k}));
        end
    end
end
