%{
TODO:
- MPC is working but constraints must be handled and circular coordinates
  must be implemented.
- Asymmetric grow
- Cartesian sampling instead of whole map

Refactor note:
This version keeps the old SNG code structure/function style as much as possible,
but replaces the graph-based logic with paper-style RSC tree growth.

Main unavoidable algorithmic changes:
- Tree is rooted at GOAL
- No overlap graph
- No Dijkstra
- Each node has exactly one parent
- Path is extracted by following parents from START-containing funnel to GOAL
%}

%%  RSC (Random Sequential Composition) With Circular Funnels

clear; clc; close all;
rng(5); % seed

% Choose scenario
scenarioId = 1;   % 1 or 2

% [map edges, obstacle polygons, start point, goal point]
[W, obs, q_start, q_goal] = getScenario(scenarioId);

%% ------------------ RSC PARAMETERS -------------------------

P.overlapThreshold = 0.1;   % kept for route intersection logic / compatibility
P.minCircArea      = 0.001;

P.maxRadius = 10.0;

P.safetyMargin  = 0.01;
P.expandStep    = 0.05;
P.enlargeStep   = 0.05;     % paper-style enlargement step
P.eta           = 0.85;     % paper suggests ~0.8 - 0.9
P.maxIter       = 5000;
P.circleRes     = 100;

%% ------------------ BUILD RSC TREE -------------------------
% polyshape, centroid, index, orientation
% parent is added for RSC tree structure
nodes = struct('poly', {}, 'c', {}, 'id', {}, 'radius', {}, 'parent', {});

% Map border polygon
workPoly = polyshape([W(1) W(2) W(2) W(1)], [W(3) W(3) W(4) W(4)]);

% Check whether the start and goal points are empty
assert(isFreePoint(q_start, obs, workPoly), 'Start is in obstacle/outside workspace.');
assert(isFreePoint(q_goal,  obs, workPoly), 'Goal is in obstacle/outside workspace.');

accepted = 0; % Number of accepted nodes

% Add START as root node
[rootPoly, rootRadius, ~] = buildCircularNode(q_start, obs, W, workPoly, P);

if isempty(rootPoly) || area(rootPoly) < P.minCircArea
    error('Root funnel at START could not be generated.');
end

accepted = accepted + 1;
nodes(accepted).poly   = rootPoly;
nodes(accepted).radius = rootRadius;
nodes(accepted).c      = q_start;
nodes(accepted).id     = accepted;
nodes(accepted).parent = 0;

startId = 1;
goalId = [];

fprintf('Added START root node as id=%d\n', startId);
fprintf('Root funnel radius = %.3f\n', nodes(startId).radius);
fprintf('Distance(start,goal) = %.3f\n', norm(q_start - q_goal));

% Trivial case: goal already inside root funnel
if isinterior(nodes(startId).poly, q_goal(1), q_goal(2))
    goalId = startId;
    fprintf('Goal is already inside the root funnel. Trivial one-funnel solution.\n');
end

iter = 0;
while iter < P.maxIter && isempty(goalId)

    iter = iter + 1;

    % Sample a free point obstacle hits are rejected
    q = sampleFreePoint(W, obs, workPoly);

    % If no point is found retry sampleFreePoint
    if isempty(q)
        continue;
    end

    % If newly sampled point is in an already covered region skip it
    if isCovered(q, nodes)
        continue;
    end

    % Find nearest existing funnel (parent)
    parentId = findNearestFunnel(q, nodes);
    qparent  = nodes(parentId).c;
    Rparent  = nodes(parentId).radius;

    % Project sample to parent boundary
    qproj = projectPointToCircleBoundary(q, qparent, Rparent);

    dir = qproj - qparent;
    nd  = norm(dir);
    if nd < 1e-12
        continue;
    end
    dir = dir / nd;

    % Paper-style new node center
    qnew = qparent + P.eta * Rparent * dir;

    if ~isFreePoint(qnew, obs, workPoly)
        continue;
    end

    % Build the funnel around the generated point
    [nodePoly0, radius0, p_contact] = buildCircularNode(qnew, obs, W, workPoly, P);

    % If created node is empty or does not meet minimum area requirement, skip
    if isempty(nodePoly0) || area(nodePoly0) < P.minCircArea
        continue;
    end

    % Enlarge funnel opposite obstacle while keeping center inside parent funnel
    [nodePoly, qnew_enl, radius] = enlargeCircularNode(qnew, radius0, p_contact, qparent, Rparent, obs, W, workPoly, P);

    if isempty(nodePoly) || area(nodePoly) < P.minCircArea
        continue;
    end

    % Accept node
    accepted = accepted + 1;
    nodes(accepted).poly   = nodePoly;
    nodes(accepted).radius = radius;
    nodes(accepted).c      = qnew_enl;
    nodes(accepted).id     = accepted;
    nodes(accepted).parent = parentId;

    fprintf('iter=%d | added node=%d | parent=%d | Rparent=%.3f | Rnew=%.3f\n', ...
        iter, accepted, parentId, Rparent, radius);

    % Stop if goal is inside newly generated funnel
    if isinterior(nodePoly, q_goal(1), q_goal(2))
        goalId = accepted;
        fprintf('Goal reached at iter=%d with %d funnels.\n', iter, accepted);
        break;
    end
end

fprintf('Terminated: accepted=%d, iter=%d (threshold=%d)\n', accepted, iter, P.maxIter);
fprintf('Accepted nodes: %d\n', accepted);

if isempty(startId)
    warning('Start was not reached within max iterations.');
end

%% ------------------ EXTRACT NODE-PATH ---------------------

% First collect [goal funnel ... start/root funnel]
pathIds = [];

if ~isempty(goalId)
    cur = goalId;
    while cur ~= 0
        pathIds(end+1) = cur; 
        cur = nodes(cur).parent;
    end
    pathIds = fliplr(pathIds);   % now [start/root funnel ... goal funnel]
end

if isempty(pathIds)
    warning('No path found in the tree.');
else
    fprintf('Found path with %d nodes.\n', numel(pathIds));
end

fprintf('startId=%d, goalId=%d\n', safeId(startId), goalId);

%% ------------------ PLOTTING -------------------------------

% RSC tree figure
figure('WindowState','maximized', 'Color','w'); hold on; axis equal;
xlim([W(1) W(2)]); ylim([W(3) W(4)]);
title('Figür 1: RSC - Tüm Kapsama Alanı (Funnels)');

% Plot obstacles
for i=1:numel(obs)
    plot(obs{i}, 'FaceColor',[0 0 0], 'FaceAlpha',0.6, 'EdgeColor','none');
end

% Plot nodes
for i=1:numel(nodes)
    plot(nodes(i).poly, 'FaceColor',[0.7 0.85 1.0], 'FaceAlpha',0.10, ...
        'EdgeColor',[0.3 0.6 1.0], 'LineWidth',0.5);

    % plot centers
    plot(nodes(i).c(1), nodes(i).c(2), '.', 'Color',[0.2 0.4 0.9], 'MarkerSize',10);

    % plot parent links
    if nodes(i).parent ~= 0
        pc = nodes(nodes(i).parent).c;
        cc = nodes(i).c;
        plot([cc(1) pc(1)], [cc(2) pc(2)], '--', 'Color',[0.5 0.5 0.5], 'LineWidth',0.7);
    end
end

% Plot start and goal points
plot(q_start(1), q_start(2), 'go', 'MarkerSize',9, 'LineWidth',2);
plot(q_goal(1),  q_goal(2),  'ro', 'MarkerSize',9, 'LineWidth',2);
grid on;

% Solution figure
figure('WindowState','maximized', 'Color','w'); hold on; axis equal;
xlim([W(1) W(2)]); ylim([W(3) W(4)]);
title('Figür 2: RSC - Seçilen Yol ve Kesişim Waypointleri');

% Plot obstacles
for i=1:numel(obs)
    plot(obs{i}, 'FaceColor',[0 0 0], 'FaceAlpha',0.6, 'EdgeColor','none');
end

route_points = [];

if ~isempty(pathIds)
    % Plot only the nodes in solution
    for k = 1:length(pathIds)
        node_idx = pathIds(k);
        plot(nodes(node_idx).poly, 'FaceColor',[1.0 0.85 0.7], ...
            'FaceAlpha',0.40, 'EdgeColor',[1.0 0.5 0.0], 'LineWidth',1.5);
    end

    num_nodes = length(pathIds);
    route_points = zeros(num_nodes + 1, 2);
    route_points(1, :) = q_start;

    for k = 1:(num_nodes - 1)
        curr_node = nodes(pathIds(k)).poly;
        next_node = nodes(pathIds(k+1)).poly;

        intersection_poly = intersect(curr_node, next_node);

        if isempty(intersection_poly) || intersection_poly.NumRegions == 0 || area(intersection_poly) <= 0
            cp = 0.5 * (nodes(pathIds(k)).c + nodes(pathIds(k+1)).c);
        else
            [cx, cy] = centroid(intersection_poly);
            cp = [cx, cy];
        end

        route_points(k+1, :) = cp;
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
    for t = 1:200
        q = [W(1) + (W(2)-W(1))*rand, W(3) + (W(4)-W(3))*rand];
        if isFreePoint(q, obs, workPoly), return; end
    end
    q = [];
end

% ok = isFreePoint(point, obstacles, map polygon)
function ok = isFreePoint(q, obs, workPoly)
    if ~isinterior(workPoly, q(1), q(2)), ok = false; return; end
    for i = 1:numel(obs)
        if isinterior(obs{i}, q(1), q(2)), ok = false; return; end
    end
    ok = true;
end

function inside = isPolyInsideWorkspace(Psh, workPoly)
    if isempty(Psh) || Psh.NumRegions == 0
        inside = false; return;
    end
    Pint = intersect(Psh, workPoly);
    inside = ~isempty(Pint) && Pint.NumRegions > 0 && abs(area(Pint) - area(Psh)) < 1e-9;
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
    for i = 1:numel(obs)
        ob = obs{i};
        if margin > 0
            try, ob = polybuffer(ob, margin); catch, end
        end
        interP = intersect(poly, ob);
        if ~isempty(interP) && interP.NumRegions > 0 && area(interP) > 0
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
            a = linspace(0,1,400)';
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

function dbox = distanceToBoxBoundary(q, W)
    dbox = min([q(1)-W(1), W(2)-q(1), q(2)-W(3), W(4)-q(2)]);
end

function [poly, radius, qobs] = buildCircularNode(qrand, obs, W, workPoly, P)
    poly = [];
    radius = 0;
    qobs = [];

    [qobs, dminObs] = closestObstaclePoint(qrand, obs);
    dminBox = distanceToBoxBoundary(qrand, W);

    % Subtract safety margin immediately so the initial poly is valid
    radius = min([dminObs, dminBox, P.maxRadius]) - P.safetyMargin;

    % Basic validity checks
    if radius <= 0 || ~isfinite(radius), return; end

    center = qrand;
    poly = circularPoly(center, radius);

    % Check if initial circle is valid
    if ~isPolyInsideWorkspace(poly, workPoly) || circCollides(poly, obs, P.safetyMargin)
        poly = [];
        radius = 0;
        return;
    end

    % Expand from the valid base
    radius = expandNode(qrand, radius, obs, W, workPoly, P);
    poly = circularPoly(qrand, radius);
end

function radius = expandNode(center, radius, obs, W, workPoly, P)
    while radius + P.expandStep <= P.maxRadius
        candR = radius + P.expandStep;

        if candR > distanceToBoxBoundary(center, W) - P.safetyMargin
            break;
        end

        cand = circularPoly(center, candR);
        if ~isPolyInsideWorkspace(cand, workPoly) || circCollides(cand, obs, P.safetyMargin)
            break;
        end

        radius = candR;
    end
end

function [poly, center, radius] = enlargeCircularNode(center0, radius0, p_contact, qparent, Rparent, obs, W, workPoly, P)
    poly   = circularPoly(center0, radius0);
    center = center0;
    radius = radius0;

    if isempty(p_contact)
        return;
    end

    v = p_contact - center0;
    nv = norm(v);
    if nv < 1e-12
        return;
    end
    v = v / nv;

    while true
        candCenter = center - P.enlargeStep * v;
        candRadius = radius + P.enlargeStep;

        % center must stay inside parent funnel
        if norm(candCenter - qparent) > Rparent
            break;
        end

        % enlarged circle must fit inside workspace
        if candRadius > distanceToBoxBoundary(candCenter, W) - P.safetyMargin
            break;
        end

        candPoly = circularPoly(candCenter, candRadius);

        if ~isPolyInsideWorkspace(candPoly, workPoly) || circCollides(candPoly, obs, P.safetyMargin)
            break;
        end

        center = candCenter;
        radius = candRadius;
        poly   = candPoly;
    end
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

function parentId = findNearestFunnel(q, nodes)
    bestVal = inf;
    parentId = 1;

    for i = 1:numel(nodes)
        val = max(norm(q - nodes(i).c) - nodes(i).radius, 0);
        if val < bestVal
            bestVal = val;
            parentId = i;
        end
    end
end

function qproj = projectPointToCircleBoundary(q, c, R)
    d = q - c;
    nd = norm(d);

    if nd < 1e-12
        qproj = c + [R 0];
    else
        qproj = c + R * d / nd;
    end
end

function out = safeId(idval)
    if isempty(idval)
        out = -1;
    else
        out = idval;
    end
end