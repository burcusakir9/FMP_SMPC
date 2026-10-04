function chain = buildChain(cfg, seed)
%BUILDCHAIN  RSC funnel chain from start to goal (RSC.m as a function).
%
%   chain = buildChain(cfg, seed)
%
% Same tree growth as RSC.m (sample a free point, skip it if covered, place
% the new center at eta*R_parent from the nearest funnel towards it, grow the
% largest circle that keeps the safety margin to obstacles and the map
% border, stop when a funnel covers the start), with the random seed as an
% input instead of rng(69). The chain is cached in cfg.dirCache.
%
% Output:
%   chain.C, chain.R   path funnel centers (2 x n) and radii (1 x n),
%                      ordered start-side -> goal (C(:,end) = q_goal)
%   chain.W, chain.obs, chain.obsAll (union of obstacles), chain.q_start,
%   chain.q_goal, chain.nNodes (tree size), chain.seed, chain.scenario

    file = fullfile(cfg.dirCache, sprintf('chain_scn%d_seed%d.mat', cfg.scenario, seed));
    if exist(file, 'file')
        s = load(file, 'chain'); chain = s.chain;
        return;
    end

    root = fileparts(fileparts(mfilename('fullpath')));
    old = cd(root);                             % getScenario reads map.kml from here
    cleanup = onCleanup(@() cd(old));
    [W, obs, q_start, q_goal] = getScenario(cfg.scenario);
    P = cfg.rsc;
    rng(seed);

    workPoly = polyshape([W(1) W(2) W(2) W(1)], [W(3) W(3) W(4) W(4)]);
    nodes = struct('poly', {}, 'c', {}, 'radius', {}, 'parent', {});

    % Goal (master) funnel
    [poly, radius] = buildCircularNode(q_goal, obs, W, workPoly, P);
    if isempty(poly) || area(poly) < P.minCircArea
        error('buildChain: goal funnel could not be generated.');
    end
    nodes(1) = struct('poly', poly, 'c', q_goal, 'radius', radius, 'parent', 0);
    startId = [];
    if isinterior(poly, q_start(1), q_start(2)), startId = 1; end

    iter = 0;
    while iter < P.maxIter && isempty(startId)
        q = sampleFreePoint(W, obs, workPoly);
        if isempty(q) || isCovered(q, nodes), continue; end

        parentId = findNearestFunnel(q, nodes);
        cp = nodes(parentId).c; Rp = nodes(parentId).radius;
        dir = q - cp;
        if norm(dir) < 1e-12, continue; end
        qnew = cp + P.eta*Rp*dir/norm(dir);
        if ~isFreePoint(qnew, obs, workPoly), continue; end

        [poly, radius] = buildCircularNode(qnew, obs, W, workPoly, P);
        if isempty(poly) || area(poly) < P.minCircArea, continue; end

        nodes(end+1) = struct('poly', poly, 'c', qnew, 'radius', radius, 'parent', parentId); %#ok<AGROW>
        if isinterior(poly, q_start(1), q_start(2))
            startId = numel(nodes);
        end
        iter = iter + 1;
    end
    if isempty(startId)
        error('buildChain: start not reached (seed %d).', seed);
    end

    % Path start -> goal
    ids = startId;
    while nodes(ids(end)).parent ~= 0
        ids(end+1) = nodes(ids(end)).parent; %#ok<AGROW>
    end

    chain.C = reshape([nodes(ids).c], 2, []);
    chain.R = [nodes(ids).radius];
    chain.W = W; chain.obs = obs; chain.obsAll = union([obs{:}]);
    chain.q_start = q_start(:); chain.q_goal = q_goal(:);
    chain.nNodes = numel(nodes); chain.seed = seed; chain.scenario = cfg.scenario;
    save(file, 'chain');
end

%% ---------------- RSC.m helper functions ----------------

function q = sampleFreePoint(W, obs, workPoly)
    for t = 1:200
        q = [W(1) + (W(2)-W(1))*rand, W(3) + (W(4)-W(3))*rand];
        if isFreePoint(q, obs, workPoly), return; end
    end
    q = [];
end

function ok = isFreePoint(q, obs, workPoly)
    ok = isinterior(workPoly, q(1), q(2));
    for i = 1:numel(obs)
        if ~ok, return; end
        ok = ~isinterior(obs{i}, q(1), q(2));
    end
end

function inside = isPolyInsideWorkspace(Psh, workPoly)
    Pint = intersect(Psh, workPoly);
    inside = Pint.NumRegions > 0 && abs(area(Pint) - area(Psh)) < 1e-9;
end

function poly = circularPoly(center, radius)
    th = linspace(0, 2*pi, 101)'; th(end) = [];
    poly = polyshape(center(1) + radius*cos(th), center(2) + radius*sin(th));
end

function coll = circCollides(poly, obs, margin)
    coll = false;
    for i = 1:numel(obs)
        ob = obs{i};
        if margin > 0
            try, ob = polybuffer(ob, margin); catch, end
        end
        I = intersect(poly, ob);
        if I.NumRegions > 0 && area(I) > 0
            coll = true; return;
        end
    end
end

function dmin = closestObstacleDistance(q, obs)
    dmin = inf;
    for i = 1:numel(obs)
        [vx, vy] = boundary(obs{i});
        if isempty(vx), continue; end
        V = [vx(:) vy(:)];
        a = linspace(0, 1, 400)';
        for k = 1:size(V,1)-1
            pts = (1-a).*V(k,:) + a.*V(k+1,:);
            dmin = min(dmin, sqrt(min(sum((pts - q).^2, 2))));
        end
    end
end

function d = distanceToBoxBoundary(q, W)
    d = min([q(1)-W(1), W(2)-q(1), q(2)-W(3), W(4)-q(2)]);
end

function [poly, radius] = buildCircularNode(q, obs, W, workPoly, P)
    poly = []; radius = 0;
    r = min([closestObstacleDistance(q, obs), distanceToBoxBoundary(q, W), P.maxRadius]) - P.safetyMargin;
    if r <= 0 || ~isfinite(r), return; end
    c = circularPoly(q, r);
    if ~isPolyInsideWorkspace(c, workPoly) || circCollides(c, obs, P.safetyMargin), return; end
    % Expand from the valid base
    while r + P.expandStep <= P.maxRadius
        rc = r + P.expandStep;
        if rc > distanceToBoxBoundary(q, W) - P.safetyMargin, break; end
        cand = circularPoly(q, rc);
        if ~isPolyInsideWorkspace(cand, workPoly) || circCollides(cand, obs, P.safetyMargin), break; end
        r = rc;
    end
    radius = r; poly = circularPoly(q, r);
end

function tf = isCovered(q, nodes)
    tf = false;
    for i = 1:numel(nodes)
        if isinterior(nodes(i).poly, q(1), q(2)), tf = true; return; end
    end
end

function id = findNearestFunnel(q, nodes)
    best = inf; id = 1;
    for i = 1:numel(nodes)
        d = abs(norm(q - nodes(i).c) - nodes(i).radius);
        if d < best, best = d; id = i; end
    end
end
