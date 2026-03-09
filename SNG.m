%{
TODO:
- MPC is working but constraints must be handles and circular coordinates
must be implemented.
- Asymmetric grow
- Cartesian sampling instead of whole map
%}

%%  SNG (Sampling-Based Neighborhood Graph) With Circular Funnels

clear; clc; close all;
rng(5); % seed

% Choose scenario
scenarioId = 1;   % 1 or 2

% [map edges, obstacle polygons, start point, goal point]
[W, obs, q_start, q_goal] = getScenario(scenarioId);

%% ------------------ SNG PARAMETERS -------------------------

P.overlapThreshold = 0.1;
P.minCircArea      = 0.001;

P.maxRadius = 10.0;

P.safetyMargin  = 0.1;
P.expandStep    = 0.05;

alpha = 0.99;
Pc    = 0.95;

m_req = ceil(log(1-Pc)/log(alpha) - 1);
m_fail = 0;                                  % consecutive failures counter

%% ------------------ BUILD SNG GRAPH ------------------------
% polyshape, centroid, index, orientation
nodes = struct('poly', {}, 'c', {}, 'id', {}, 'radius', {});

% Map border polygon
workPoly = polyshape([W(1) W(2) W(2) W(1)],[W(3) W(3) W(4) W(4)]);

% Check wheter the start and goal points are empty
assert(isFreePoint(q_start, obs, workPoly), 'Start is in obstacle/outside workspace.');
assert(isFreePoint(q_goal,  obs, workPoly), 'Goal is in obstacle/outside workspace.');

accepted = 0; % Number of accepted nodes

m = 0;     
while m < m_req

    % Sample a free point obstacle hits are rejected
    q = sampleFreePoint(W, obs, workPoly);

    % If no point is found retry sampleFreePoint 
    if isempty(q)
        continue; 
    end

    % If newly sampled point is in a already covered region increase failure counter
    if isCovered(q, nodes)
        m = m + 1;
        continue;
    end

    % Build the funnel around the sampled point
    [nodePoly, radius] = buildCircularNode(q, obs, workPoly, P);

    % If created node is empty or does not meet minimum area requirement, skip without accepting
    if isempty(nodePoly) || area(nodePoly) < P.minCircArea
        continue;
    end

    % Accept node 
    accepted = accepted + 1;
    nodes(accepted).poly  = nodePoly;
    nodes(accepted).radius = radius;
    [cx, cy] = centroid(nodePoly);
    nodes(accepted).c  = [cx, cy];
    nodes(accepted).id = accepted;

    % Reset failure counter after success
    m = 0;
end

fprintf('Terminated: accepted=%d, m=%d (threshold=%d)\n', accepted, m, m_req);
fprintf('Accepted nodes: %d\n', accepted);

% Add start and goal nodes 
nodes = addPointAsNode(nodes, q_start, obs, workPoly, P, "START");
startId = numel(nodes);

nodes = addPointAsNode(nodes, q_goal, obs, workPoly, P, "GOAL");
goalId = numel(nodes);

% A(i,j) = weight if node i and j overlap
% A(i,j) = 0 if they don’t
A = rebuildAdjacency(nodes, P);

%% ------------------ SHORTEST NODE-PATH ---------------------

[pathIds, distVal] = dijkstraSparse(A, goalId, startId); %% Search fromn goal to start

pathIds = flip(pathIds); %% Flip the waypoints.

if isempty(pathIds)
    warning('No path found in the graph (graph may be disconnected).');
else
    fprintf('Found path with %d nodes, total centroid-distance = %.3f\n', numel(pathIds), distVal);
end

fprintf('deg(start)=%d, deg(goal)=%d\n', nnz(A(startId,:)), nnz(A(goalId,:)));

%% ------------------ PLOTTING -------------------------------

% SNG figure
figure('WindowState','maximized', 'Color','w'); hold on; axis equal;
xlim([W(1) W(2)]); ylim([W(3) W(4)]);
title('Figür 1: SNG - Tüm Kapsama Alanı (Funnels)');

% Plot obstacles
for i=1:numel(obs)
    plot(obs{i}, 'FaceColor',[0 0 0], 'FaceAlpha',0.6, 'EdgeColor','none');
end

% Plot nodes
for i=1:numel(nodes)
    plot(nodes(i).poly, 'FaceColor',[0.7 0.85 1.0], 'FaceAlpha',0.10, 'EdgeColor',[0.3 0.6 1.0], 'LineWidth',0.5);
end

% Plot start and goal points
plot(q_start(1), q_start(2), 'go', 'MarkerSize',9, 'LineWidth',2);
plot(q_goal(1),  q_goal(2),  'ro', 'MarkerSize',9, 'LineWidth',2);
grid on;

% Dijkstra figure
figure('WindowState','maximized', 'Color','w'); hold on; axis equal;
xlim([W(1) W(2)]); ylim([W(3) W(4)]);
title('Figür 2: Dijkstra - Seçilen Yol ve Kesişim Waypointleri');

% Plot obstacles
for i=1:numel(obs)
    plot(obs{i}, 'FaceColor',[0 0 0], 'FaceAlpha',0.6, 'EdgeColor','none');
end

if ~isempty(pathIds)
    % Plot only the nodes in Dijkstra solutionn
    for k = 1:length(pathIds)
        node_idx = pathIds(k);
        plot(nodes(node_idx).poly, 'FaceColor',[1.0 0.85 0.7], 'FaceAlpha',0.40, 'EdgeColor',[1.0 0.5 0.0], 'LineWidth',1.5);
    end
    
    num_nodes = length(pathIds);
    route_points = zeros(num_nodes + 1, 2); 
    route_points(1, :) = q_start; 
    
    for k = 1:(num_nodes - 1)
        curr_node = nodes(pathIds(k)).poly;
        next_node = nodes(pathIds(k+1)).poly;
        
        intersection_poly = intersect(curr_node, next_node);
        [cx, cy] = centroid(intersection_poly);
        route_points(k+1, :) = [cx, cy];
    end
    
    route_points(end, :) = q_goal; 
    
    % Plot cross points
    plot(route_points(:,1), route_points(:,2), 'bo', 'MarkerSize',6, 'LineWidth',1.5);
end

% Plot start and goal points
plot(q_start(1), q_start(2), 'go', 'MarkerSize',9, 'LineWidth',2);
plot(q_goal(1),  q_goal(2),  'ro', 'MarkerSize',9, 'LineWidth',2);
grid on;
%% ===================== FUNCTIONS ===========================

% point = sampleFreePoint(map corner points, obstacles, map polygon)
function q = sampleFreePoint(W, obs, workPoly)
    for t=1:200
        q = [W(1) + (W(2)-W(1))*rand, W(3) + (W(4)-W(3))*rand];
        if isFreePoint(q, obs, workPoly), return; end
    end
    q = [];
end

% ok = isFreePoint(point, obstacles, map polygon)
function ok = isFreePoint(q, obs, workPoly)
    if ~isinterior(workPoly, q(1), q(2)), ok=false; return; end
    for i=1:numel(obs)
        if isinterior(obs{i}, q(1), q(2)), ok=false; return; end
    end
    ok = true;
end

function inside = isPolyInsideWorkspace(Psh, workPoly)
    if isempty(Psh) || Psh.NumRegions == 0
        inside = false; return;
    end
    Pint = intersect(Psh, workPoly);
    inside = ~isempty(Pint) && abs(area(Pint) - area(Psh)) < 1e-9;
end


function poly = circularPoly(center, radius)
    numPoints = 100; % Resolution of the circle
    theta = linspace(0, 2*pi, numPoints)';
    x = center(1) + radius * cos(theta);
    y = center(2) + radius * sin(theta);
    poly = polyshape(x, y);
end

function coll = circCollides(poly, obs, margin)
    coll = false;
    for i=1:numel(obs)
        ob = obs{i};
        if margin > 0
            try, ob = polybuffer(ob, margin); catch, end
        end
        interP = intersect(poly, ob);
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

function [poly, radius] = buildCircularNode(qrand, obs, workPoly, P)
    poly = [];
    radius = 0;
    [~, dmin] = closestObstaclePoint(qrand, obs);
    
    % Subtract safety margin immediately so the initial poly is valid
    radius = dmin - P.safetyMargin;
    
    % Basic validity checks
    if radius <= 0 || ~isfinite(radius), return; end
    if radius > P.maxRadius, radius = P.maxRadius; end
    
    center = qrand;
    poly = circularPoly(center, radius);
    
    % Check if initial circle is valid
    if ~isPolyInsideWorkspace(poly, workPoly) || circCollides(poly, obs, P.safetyMargin)
        poly = []; 
        return; 
    end
    
    % Expand from the valid base
    radius = expandNode(qrand, radius, obs, workPoly, P);
    poly = circularPoly(qrand, radius);
end


function radius = expandNode(center, radius, obs, workPoly, P)
    while radius + P.expandStep <= P.maxRadius
        cand = circularPoly(center, radius);
        if ~isPolyInsideWorkspace(cand, workPoly) || circCollides(cand, obs, P.safetyMargin), break; end
        radius = radius + P.expandStep;
    end

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
    [poly, radius] = buildCircularNode(q, obs, workPoly, P);
    
    % Fallback: if build failed, try a tiny circle
    if isempty(poly)
        radius = P.expandStep; 
        poly = circularPoly(q, radius);
    end
    
    % Now check validity
    if ~isPolyInsideWorkspace(poly, workPoly) % || circCollides(poly, obs, P.safetyMargin)
         error('%s node could not be embedded. Point is likely too close to obstacle.', tag);
    end

    % ensure START/GOAL overlaps at least one existing node
    if ~isempty(nodes)
        if ~hasAnyOverlap(poly, nodes, P.overlapThreshold)
            % grow uniformly a bit until it overlaps something (still collision-free)
            [poly, ok] = growUntilOverlap(poly, q, radius, nodes, obs, workPoly, P);
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
    nodes(n).radius = radius;

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

function [poly, ok] = growUntilOverlap(poly, center, radius, nodes, obs, workPoly, P)
    ok = false;

    for k=1:200
        radius = radius + P.expandStep;
        cand = circularPoly(center, radius);
        if ~isPolyInsideWorkspace(cand, workPoly) || circCollides(cand, obs, P.safetyMargin)
            break;
        end
        poly = cand;
        if hasAnyOverlap(poly, nodes, P.overlapThreshold)
            ok = true; return;
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

% Check if the point is in any of the nodes
function tf = isCovered(q, nodes)
    tf = false;
    for i = 1:numel(nodes)
        if isinterior(nodes(i).poly, q(1), q(2))
            tf = true;
            return;
        end
    end
end
