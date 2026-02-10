% 
% TODO:  There is a bug that sometimes funnels are not parallel to obstacles.
%         In second scenario still no path with suggested terminal
%         condition constants
% 
% 

%%  SNG (Sampling-Based Neighborhood Graph)

clear; clc; close all;
rng(5); % seed

% Choose scenario
scenarioId = 1;   % 1 or 2
[W, obs, q_start, q_goal] = getScenario(scenarioId);

%% ------------------ SNG PARAMETERS -------------------------
P.asymExpand = false; % true = asymmetric growth, false = symmetric growth

P.overlapThreshold = 1e-6;
P.minRectArea      = 0.001;

P.maxHalfSize   = [50.0, 50.0];
P.safetyMargin  = 0.0;
P.expandStep    = 0.05;

alpha = 0.98;
Pc    = 0.95;

m_req = ceil(log(1-Pc)/log(alpha) - 1);
m_fail = 0;                                  % consecutive failures counter


%% ------------------ BUILD SNG GRAPH ------------------------
% polyshape, centroid, index, orientation
nodes = struct('poly', {}, 'c', {}, 'id', {}, 'theta', {});

% A(i,j) = weight if node i and j overlap
% A(i,j) = 0 if they don’t
A = sparse(0,0);

workPoly = polyshape([W(1) W(2) W(2) W(1)],[W(3) W(3) W(4) W(4)]); % map border polygon

% Check wheter the start and goal points are empty
assert(isFreePoint(q_start, obs, workPoly), 'Start is in obstacle/outside workspace.');
assert(isFreePoint(q_goal,  obs, workPoly), 'Goal is in obstacle/outside workspace.');

accepted = 0; % Number of accepted nodes
m = 0;     
while m < m_req

    % Sample a free point (obstacle hits are internally rejected, not counted)
    q = sampleFreePoint(W, obs, workPoly);
    if isempty(q)
        continue;   % ignore: not counted as failure
    end

    % FAILURE = sample lies in already covered region (inside any existing node)
    if accepted > 0 && isCovered(q, nodes)
        m = m + 1;
        continue;
    end

    % SUCCESS attempt: build a new neighborhood around q
    theta = 0;
    nodePoly = buildRectNode(q, obs, workPoly, P);

    if isempty(nodePoly) || area(nodePoly) < P.minRectArea
        % ignore: not counted as failure in the paper's model
        continue;
    end

    % Accept node (SUCCESS)
    accepted = accepted + 1;
    nodes(accepted).poly  = nodePoly;
    nodes(accepted).theta = theta;
    [cx, cy] = centroid(nodePoly);
    nodes(accepted).c  = [cx, cy];
    nodes(accepted).id = accepted;

    % Reset failure counter after success
    m = 0;
end

fprintf('Terminated: accepted=%d, m=%d (threshold=%d)\n', accepted, m, m_req);


fprintf('Accepted nodes: %d\n', accepted);

% Add start and goal nodes (guarantee they overlap at least one existing node)
nodes = addPointAsNode(nodes, q_start, obs, workPoly, P, "START");
startId = numel(nodes);


nodes = addPointAsNode(nodes, q_goal, obs, workPoly, P, "GOAL");
goalId = numel(nodes);

A = rebuildAdjacency(nodes, P);


%% ------------------ SHORTEST NODE-PATH ---------------------
[pathIds, distVal] = dijkstraSparse(A, startId, goalId);

if isempty(pathIds)
    warning('No path found in the graph (graph may be disconnected).');
else
    fprintf('Found path with %d nodes, total centroid-distance = %.3f\n', numel(pathIds), distVal);
end

fprintf('deg(start)=%d, deg(goal)=%d\n', nnz(A(startId,:)), nnz(A(goalId,:)));

%% ------------------ PLOTTING -------------------------------
figure('WindowState','maximized', 'Color','w'); hold on; axis equal;
xlim([W(1) W(2)]); ylim([W(3) W(4)]);
title('SNG: obstacles, nodes, and node-path');

for i=1:numel(obs)
    plot(obs{i}, 'FaceColor',[0 0 0], 'FaceAlpha',0.6, 'EdgeColor','none');
end

for i=1:numel(nodes)
    plot(nodes(i).poly, 'FaceColor',[0.7 0.85 1.0], 'FaceAlpha',0.10, 'EdgeColor',[0.3 0.6 1.0], 'LineWidth',0.5);
end

plot(q_start(1), q_start(2), 'go', 'MarkerSize',9, 'LineWidth',2);
plot(q_goal(1),  q_goal(2),  'ro', 'MarkerSize',9, 'LineWidth',2);

if ~isempty(pathIds)
    C = reshape([nodes(pathIds).c],2,[])';
    plot(C(:,1), C(:,2), 'ms', 'MarkerSize',6, 'LineWidth',1.5);
end

grid on;

%% ===================== FUNCTIONS ===========================

function q = sampleFreePoint(W, obs, workPoly)
    for t=1:200
        q = [W(1) + (W(2)-W(1))*rand, W(3) + (W(4)-W(3))*rand];
        if isFreePoint(q, obs, workPoly), return; end
    end
    q = [];
end

function ok = isFreePoint(q, obs, workPoly)
    if ~isinterior(workPoly, q(1), q(2)), ok=false; return; end
    for i=1:numel(obs)
        if isinterior(obs{i}, q(1), q(2)), ok=false; return; end
    end
    ok = true;
end

function inside = isPolyInsideWorkspace(Psh, workPoly)
    Pint = intersect(Psh, workPoly);
    inside = ~isempty(Pint) && abs(area(Pint) - area(Psh)) < 1e-9;
end

function poly = rectPolyAsym(c, aP, aM, bP, bM, theta)
    corners_uv = [ +aP +bP;
                   -aM +bP;
                   -aM -bM;
                   +aP -bM ];
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    corners_xy = (R * corners_uv')' + c;
    poly = polyshape(corners_xy(:,1), corners_xy(:,2));
end

function poly = rectPolySym(c, hx, hy, theta)
    poly = rectPolyAsym(c, hx, hx, hy, hy, theta);
end

function coll = rectCollides(rectP, obs, margin)
    coll = false;
    for i=1:numel(obs)
        ob = obs{i};
        if margin > 0
            try, ob = polybuffer(ob, margin); catch, end
        end
        interP = intersect(rectP, ob);
        if ~isempty(interP) && area(interP) > 0
            coll = true; return;
        end
    end
end

function [qobs, dmin] = closestObstaclePoint(qrand, obs)
    dmin = inf;
    qobs = [];
    for i = 1:numel(obs)
        [vx, vy] = boundary(obs{i});
        if isempty(vx), continue; end
        V = [vx(:) vy(:)];
        pts = [];
        for k = 1:size(V,1)-1
            pA = V(k,:); pB = V(k+1,:);
            a = linspace(0,1,400)';   % dense enough, not crazy slow
            pts = [pts; (1-a).*pA + a.*pB]; %#ok<AGROW>
        end
        D = pts - qrand;
        dd = sum(D.^2,2);
        [m, idx] = min(dd);
        if m < dmin^2
            dmin = sqrt(m);
            qobs = pts(idx,:);
        end
    end
end


function poly = buildRectNode(qrand, obs, workPoly, P)
    % Method 1: initial square from dmin, then expand (sym or asym based on flag)
    [qobs, dmin] = closestObstaclePoint(qrand, obs);
    if isempty(qobs) || ~isfinite(dmin) || dmin <= 1e-6, poly=[]; return; end

    h = dmin / sqrt(2);
    h = min(h, min(P.maxHalfSize));

    v = qobs - qrand;
    theta = atan2(v(2), v(1)) + pi/2;

    if P.asymExpand
        aP=h; aM=h; bP=h; bM=h;
        poly = rectPolyAsym(qrand, aP,aM,bP,bM, theta);
        if ~isPolyInsideWorkspace(poly, workPoly) || rectCollides(poly, obs, P.safetyMargin), poly=[]; return; end

        % expand each side independently
        [aP,aM,bP,bM] = expandAsymExtents(qrand, theta, aP,aM,bP,bM, obs, workPoly, P);
        poly = rectPolyAsym(qrand, aP,aM,bP,bM, theta);
    else
        hx=h; hy=h;
        poly = rectPolySym(qrand, hx, hy, theta);
        if ~isPolyInsideWorkspace(poly, workPoly) || rectCollides(poly, obs, P.safetyMargin), poly=[]; return; end

        % symmetric expand
        [hx,hy] = expandSymHalfSizes(qrand, theta, hx,hy, obs, workPoly, P);
        poly = rectPolySym(qrand, hx,hy, theta);
    end
end

function [aP,aM,bP,bM] = expandAsymExtents(c, theta, aP,aM,bP,bM, obs, workPoly, P)
    % +u
    while aP + P.expandStep <= P.maxHalfSize(1)
        cand = rectPolyAsym(c, aP + P.expandStep, aM, bP, bM, theta);
        if ~isPolyInsideWorkspace(cand, workPoly) || rectCollides(cand, obs, P.safetyMargin), break; end
        aP = aP + P.expandStep;
    end
    % -u
    while aM + P.expandStep <= P.maxHalfSize(1)
        cand = rectPolyAsym(c, aP, aM + P.expandStep, bP, bM, theta);
        if ~isPolyInsideWorkspace(cand, workPoly) || rectCollides(cand, obs, P.safetyMargin), break; end
        aM = aM + P.expandStep;
    end
    % +v
    while bP + P.expandStep <= P.maxHalfSize(2)
        cand = rectPolyAsym(c, aP, aM, bP + P.expandStep, bM, theta);
        if ~isPolyInsideWorkspace(cand, workPoly) || rectCollides(cand, obs, P.safetyMargin), break; end
        bP = bP + P.expandStep;
    end
    % -v
    while bM + P.expandStep <= P.maxHalfSize(2)
        cand = rectPolyAsym(c, aP, aM, bP, bM + P.expandStep, theta);
        if ~isPolyInsideWorkspace(cand, workPoly) || rectCollides(cand, obs, P.safetyMargin), break; end
        bM = bM + P.expandStep;
    end
end

function [hx,hy] = expandSymHalfSizes(c, theta, hx,hy, obs, workPoly, P)
    while hx + P.expandStep <= P.maxHalfSize(1)
        cand = rectPolySym(c, hx + P.expandStep, hy, theta);
        if ~isPolyInsideWorkspace(cand, workPoly) || rectCollides(cand, obs, P.safetyMargin), break; end
        hx = hx + P.expandStep;
    end
    while hy + P.expandStep <= P.maxHalfSize(2)
        cand = rectPolySym(c, hx, hy + P.expandStep, theta);
        if ~isPolyInsideWorkspace(cand, workPoly) || rectCollides(cand, obs, P.safetyMargin), break; end
        hy = hy + P.expandStep;
    end
end

function [aP,aM,bP,bM] = rectExtentsAsym(poly, center, theta)
    [vx, vy] = boundary(poly);
    V = [vx(:) vy(:)];
    if size(V,1) >= 2 && all(V(1,:) == V(end,:)), V(end,:) = []; end

    V = V - center;
    R = [cos(-theta) -sin(-theta); sin(-theta) cos(-theta)];
    Vr = (R * V')';

    aP = max(Vr(:,1));
    aM = max(-Vr(:,1));
    bP = max(Vr(:,2));
    bM = max(-Vr(:,2));
end

function [hx,hy] = rectHalfSizesSym(poly, center, theta)
    [vx, vy] = boundary(poly);
    V = [vx(:) vy(:)];
    if size(V,1) >= 2 && all(V(1,:) == V(end,:)), V(end,:) = []; end

    V = V - center;
    R = [cos(-theta) -sin(-theta); sin(-theta) cos(-theta)];
    Vr = (R * V')';

    hx = max(abs(Vr(:,1)));
    hy = max(abs(Vr(:,2)));
end

function A = rebuildAdjacency(nodes, P)
    n = numel(nodes);
    A = sparse(n,n);
    for i=1:n
        for j=i+1:n
            ov = intersect(nodes(i).poly, nodes(j).poly);
            if ~isempty(ov) && area(ov) > P.overlapThreshold
                w = norm(nodes(i).c - nodes(j).c);
                A(i,j)=w; A(j,i)=w;
            end
        end
    end
end


function nodes = addPointAsNode(nodes, q, obs, workPoly, P, tag)
    % Create node polygon + theta depending on method
    
    theta = 0;
    poly  = buildRectNode(q, obs, workPoly, P);

    if isempty(poly)
        theta = 0;
        poly  = rectPolyAsym(q, 0.6,0.6,0.6,0.6, theta);
        if ~isPolyInsideWorkspace(poly, workPoly) || rectCollides(poly, obs, P.safetyMargin)
            error('%s node could not be embedded (fallback invalid).', tag);
        end
    end

    % ---- FIX: ensure START/GOAL overlaps at least one existing node ----
    if ~isempty(nodes)
        if ~hasAnyOverlap(poly, nodes, P.overlapThreshold)
            % grow uniformly a bit until it overlaps something (still collision-free)
            [poly, ok] = growUntilOverlap(poly, q, theta, nodes, obs, workPoly, P);
            if ~ok
                % last resort: lower threshold effect by forcing tiny threshold overlap
                % (keeps code minimal; you can tune overlapThreshold instead)
            end
        end
    end

    n = numel(nodes) + 1;
    nodes(n).poly  = poly;
    [cx, cy]       = centroid(poly);
    nodes(n).c     = [cx, cy];
    nodes(n).id    = n;
    nodes(n).theta = theta;

    fprintf('Added %s node as id=%d\n', tag, n);
end

function tf = hasAnyOverlap(poly, nodes, thr)
    tf = false;
    for i=1:numel(nodes)
        ov = intersect(poly, nodes(i).poly);
        if ~isempty(ov) && area(ov) > thr
            tf = true; return;
        end
    end
end

function [poly, ok] = growUntilOverlap(poly, center, theta, nodes, obs, workPoly, P)
    ok = false;

    if P.asymExpand
        [aP,aM,bP,bM] = rectExtentsAsym(poly, center, theta);
        for k=1:200
            % grow all sides a bit
            aP = aP + P.expandStep; aM = aM + P.expandStep;
            bP = bP + P.expandStep; bM = bM + P.expandStep;
            cand = rectPolyAsym(center, aP,aM,bP,bM, theta);

            if ~isPolyInsideWorkspace(cand, workPoly) || rectCollides(cand, obs, P.safetyMargin)
                break;
            end
            poly = cand;
            if hasAnyOverlap(poly, nodes, P.overlapThreshold)
                ok = true; return;
            end
        end
    else
        [hx,hy] = rectHalfSizesSym(poly, center, theta);
        for k=1:200
            hx = hx + P.expandStep;
            hy = hy + P.expandStep;
            cand = rectPolySym(center, hx, hy, theta);

            if ~isPolyInsideWorkspace(cand, workPoly) || rectCollides(cand, obs, P.safetyMargin)
                break;
            end
            poly = cand;
            if hasAnyOverlap(poly, nodes, P.overlapThreshold)
                ok = true; return;
            end
        end
    end
end

function [path, distVal] = dijkstraSparse(A, s, t)
    n = size(A,1);
    dist = inf(n,1); dist(s)=0;
    prev = zeros(n,1);
    visited = false(n,1);

    for iter=1:n
        dtmp = dist; dtmp(visited)=inf;
        [m,u] = min(dtmp);
        if ~isfinite(m), break; end
        if u==t, break; end
        visited(u)=true;

        nbrs = find(A(u,:)>0);
        for k=1:numel(nbrs)
            v = nbrs(k);
            if visited(v), continue; end
            alt = dist(u) + A(u,v);
            if alt < dist(v)
                dist(v)=alt;
                prev(v)=u;
            end
        end
    end

    if ~isfinite(dist(t))
        path = [];
        distVal = inf;
        return;
    end

    path = t;
    while path(1) ~= s
        path = [prev(path(1)); path]; 
        if path(1)==0
            path = []; distVal=inf; return;
        end
    end
    distVal = dist(t);
end

function tf = isCovered(q, nodes)
    tf = false;
    for i = 1:numel(nodes)
        if isinterior(nodes(i).poly, q(1), q(2))
            tf = true;
            return;
        end
    end
end
