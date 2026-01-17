%%  SNG (Sampling-Based Neighborhood Graph

clear; 
clc; 
close all;
rng(4);  % random number generator

%% ------------------ MAP DEFINITION -------------------------
% Workspace bounds 
W = [0 30 0 20];  % [xmin xmax ymin ymax]

% Obstacles as polynoms
obs = {};
obs{end+1} = polyshape([6 10 10 6],[3 3 7 7]); % [x1 x2 x3 x4],[y1 y2 y3 y4]
obs{end+1} = polyshape([14 18 18 14],[2 2 6 6]);
obs{end+1} = polyshape([20 24 24 20],[10 10 16 16]);
obs{end+1} = polyshape([8 12 12 8],[12 12 18 18]);
obs{end+1} = polyshape([2 4 4 2],[9 9 14 14]);

% Start / Goal points [x, y]
q_start = [2 2]; 
q_goal  = [28 18];

%% ------------------ SNG PARAMETERS -------------------------
P.maxNodes            = 100;    % number of nodes to attempt/accept
P.maxTriesPerNode     = 10;     % attempts before giving up sampling
P.overlapThreshold    = 1e-3;   % overlap area threshold for edge creation
P.minRectArea         = 0.40;   % reject tiny nodes
P.seedRectHalfSize    = [0.6 0.6]; % initial half-widths [hx hy]
P.maxExpandSteps      = 18;     % rectangle expansion iterations
P.expandGrowth        = 1.25;   % multiply half-sizes by this if safe
P.maxHalfSize         = [6 6];  % cap rectangle half-sizes
P.safetyMargin        = 0.05;   % Safety margin so rectangles don't graze obstacles too tightly

%% ------------------ BUILD SNG GRAPH ------------------------
% Node struct:
% node(i).poly   polyshape
% node(i).c      centroid [x y]
% node(i).id     integer index

nodes = struct('poly', {}, 'c', {}, 'id', {});
A = sparse(0,0); % adjacency (weighted)

% Workspace polygon
workPoly = polyshape([W(1) W(2) W(2) W(1)],[W(3) W(3) W(4) W(4)]);

% Validate start/goal
assert(isFreePoint(q_start, obs, workPoly), 'Start is in obstacle/outside workspace.');
assert(isFreePoint(q_goal,  obs, workPoly), 'Goal is in obstacle/outside workspace.');

% Build random nodes
accepted = 0;

while accepted < P.maxNodes
    success = false;

    for t = 1:P.maxTriesPerNode
        q = sampleFreePoint(W, obs, workPoly);
        if isempty(q), continue; end

        nodePoly = buildRectNode(q, obs, workPoly, P);
        if isempty(nodePoly), continue; end

        % Check is minimum area threshold holds
        if area(nodePoly) < P.minRectArea
            continue;
        end

        % ---- accept node ----

        % iterate
        accepted = accepted + 1;

        % Add accepted node to the list
        nodes(accepted).poly = nodePoly;

        % Obtain node center
        [cx, cy] = centroid(nodePoly);

        % Record node centroid info
        nodes(accepted).c = [cx, cy];

        % Record node id
        nodes(accepted).id = accepted;

        success = true;
        break;
    end

    % If couldn't create a new node in maxTriesPerNode attempts, stop
    if ~success
        break;
    end
end

fprintf('Accepted nodes: %d\n', accepted);

% Add start and goal as nodes and connect by overlap
nodes = addPointAsNode(nodes, q_start, obs, workPoly, P, "START");
nodes = addPointAsNode(nodes, q_goal,  obs, workPoly, P, "GOAL");

% Rebuild adjacency including the new nodes
A = rebuildAdjacency(nodes, P);

startId = findNodeContainingPoint(nodes, q_start);
goalId  = findNodeContainingPoint(nodes, q_goal);

if isempty(startId) || isempty(goalId)
    error('Start/Goal could not be embedded into the graph. Try increasing maxNodes or changing parameters.');
end

%% ------------------ SHORTEST NODE-PATH ---------------------

% Find the shortest path
[pathIds, distVal] = dijkstraSparse(A, startId, goalId);  % TODO start goal points are not in the path

if isempty(pathIds)
    warning('No path found in the graph (graph may be disconnected).');
else
    fprintf('Found path with %d nodes, total centroid-distance = %.3f\n', numel(pathIds), distVal);
end

%% ------------------ PLOTTING -------------------------------
figure('WindowState','maximized', 'Color','w'); hold on; axis equal;
xlim([W(1) W(2)]); ylim([W(3) W(4)]);
title('SNG: obstacles, rectangular nodes, overlap edges, and node-path');

% Workspace boundary
% plot(workPoly, 'FaceColor','none','EdgeColor',[0 0 0],'LineWidth',1.2);

% Obstacles
for i=1:numel(obs)
    plot(obs{i}, 'FaceColor',[0.0 0.0 0.0], 'FaceAlpha',0.6, 'EdgeColor','none');
end

% Nodes (light)
for i=1:numel(nodes)
    plot(nodes(i).poly, 'FaceColor',[0.7 0.85 1.0], 'FaceAlpha',0.10, 'EdgeColor',[0.3 0.6 1.0], 'LineWidth',0.5);
end

% Graph edges
[idx_i, idx_j, w] = find(triu(A,1));
for k=1:numel(w)
    ci = nodes(idx_i(k)).c;
    cj = nodes(idx_j(k)).c;
    % plot([ci(1) cj(1)], [ci(2) cj(2)], '-', 'Color',[0.6 0.6 0.6], 'LineWidth',0.6);
end

% Start/Goal points
plot(q_start(1), q_start(2), 'go', 'MarkerSize',9, 'LineWidth',2);
plot(q_goal(1),  q_goal(2),  'ro', 'MarkerSize',9, 'LineWidth',2);

% Path (draw centroid-to-centroid)
if ~isempty(pathIds)
    C = reshape([nodes(pathIds).c],2,[])';
    plot(C(:,1), C(:,2), 'm-', 'LineWidth',3);
    plot(C(:,1), C(:,2), 'ms', 'MarkerSize',6, 'LineWidth',1.5);
end

legend({'Obstacles','Nodes','Start','Goal','Node-path'}, 'Location','northeastoutside');
grid on;

%% ===================== FUNCTIONS ===========================

function q = sampleFreePoint(W, obs, workPoly)
    % Uniform sample in workspace bounds, reject if in obstacle/outside.
    for t=1:200
        q = [W(1) + (W(2)-W(1))*rand, W(3) + (W(4)-W(3))*rand];
        if isFreePoint(q, obs, workPoly)
            return;
        end
    end
    q = [];
end

function ok = isFreePoint(q, obs, workPoly)
    if ~isinterior(workPoly, q(1), q(2))
        ok = false; return;
    end
    for i=1:numel(obs)
        if isinterior(obs{i}, q(1), q(2))
            ok = false; return;
        end
    end
    ok = true;
end

function poly = buildRectNode(q, obs, workPoly, P)
    theta = 0;
    [dirVec, dmin] = nearestObstacleDirection(q, obs);
    if ~isempty(dirVec) && isfinite(dmin) && dmin > 1e-6
        theta = atan2(dirVec(2), dirVec(1)) + pi/2;
    end


    hx = P.seedRectHalfSize(1);
    hy = P.seedRectHalfSize(2);

    % Initial candidate (NO intersect/clip)
    poly0 = rectPoly(q, hx, hy, theta);

    % Reject if outside workspace or colliding
    if ~isPolyInsideWorkspace(poly0, workPoly) || rectCollides(poly0, obs, P.safetyMargin)
        poly = []; return;
    end

    poly = poly0;

    for k = 1:P.maxExpandSteps
        hx2 = min(hx * P.expandGrowth, P.maxHalfSize(1));
        hy2 = min(hy * P.expandGrowth, P.maxHalfSize(2));

        cand = rectPoly(q, hx2, hy2, theta);

        % Stop expansion if it would leave workspace
        if ~isPolyInsideWorkspace(cand, workPoly)
            break;
        end

        if ~rectCollides(cand, obs, P.safetyMargin)
            poly = cand; hx = hx2; hy = hy2;
        else
            % Try anisotropic expansion
            cand1 = rectPoly(q, hx2, hy, theta);
            cand2 = rectPoly(q, hx,  hy2, theta);

            ok1 = isPolyInsideWorkspace(cand1, workPoly) && ~rectCollides(cand1, obs, P.safetyMargin);
            ok2 = isPolyInsideWorkspace(cand2, workPoly) && ~rectCollides(cand2, obs, P.safetyMargin);

            if ok1 && area(cand1) >= area(poly)
                poly = cand1; hx = hx2;
            end
            if ok2 && area(cand2) >= area(poly)
                poly = cand2; hy = hy2;
            end

            if ~(ok1 || ok2)
                break;
            end
        end
    end
end


function inside = isPolyInsideWorkspace(P, workPoly)
    Pint = intersect(P, workPoly);
    inside = ~isempty(Pint) && abs(area(Pint) - area(P)) < 1e-9;
end



function poly = rectPoly(c, hx, hy, theta)
    % Rectangle centered at c with half-sizes hx, hy, rotated by theta.
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    corners = [ -hx -hy;
                 hx -hy;
                 hx  hy;
                -hx  hy ];
    corners = (R * corners')';
    corners = corners + c;
    poly = polyshape(corners(:,1), corners(:,2));
end

function coll = rectCollides(rectP, obs, margin)
    % Collision if intersects any obstacle (with small margin)
    % Margin implemented by slightly buffering obstacles outward (approx via polybuffer if available).
    coll = false;
    for i=1:numel(obs)
        ob = obs{i};
        if margin > 0
            try
                ob = polybuffer(ob, margin);
            catch
                % if polybuffer not available, ignore margin
            end
        end
        interP = intersect(rectP, ob);
        if ~isempty(interP) && area(interP) > 0
            coll = true; return;
        end
    end
end

function [dirVec, dmin] = nearestObstacleDirection(q, obs)
    % Approximate nearest obstacle direction by scanning obstacle vertices.
    % Returns vector from q toward nearest obstacle vertex.
    dmin = inf;
    dirVec = [];
    for i=1:numel(obs)
        [vx, vy] = boundary(obs{i});
        if isempty(vx); continue; end
        V = [vx(:) vy(:)];
        D = V - q;
        dd = sum(D.^2,2);
        [m, idx] = min(dd);
        if m < dmin
            dmin = m;
            dirVec = D(idx,:);
        end
    end
    dmin = sqrt(dmin);
    if isempty(dirVec), dirVec=[]; end
end

function nodes = addPointAsNode(nodes, q, obs, workPoly, P, tag)
    % Add a node at q so start/goal can attach by overlap.
    P2 = P;
    P2.seedRectHalfSize = [0.8 0.8];  % TODO this shoud not be fixed size
    P2.maxExpandSteps   = 6;
    poly = buildRectNode(q, obs, workPoly, P2);

    if isempty(poly)
    % fallback: tiny unrotated rect (NO clipping)
    poly = rectPoly(q, 0.6, 0.6, 0);
        if ~isPolyInsideWorkspace(poly, workPoly) || rectCollides(poly, obs, P.safetyMargin)
            error('%s node could not be embedded (fallback invalid).', tag);
        end
    end

    n = numel(nodes) + 1;
    nodes(n).poly = poly;
    [cx, cy] = centroid(poly);
    nodes(n).c = [cx, cy];
    nodes(n).id = n;
    fprintf('Added %s node as id=%d\n', tag, n);
end

function A = rebuildAdjacency(nodes, P)
    n = numel(nodes); % # of nodes
    A = sparse(n,n);  % Square matrix by the size of # of nodes

    % Iterate over all node pairs
    for i=1:n
        for j=i+1:n

            % Centroid of the node pairs
            ci = nodes(i).c; 
            cj = nodes(j).c;
    
            % Overlap polygon of node pairs
            ov = intersect(nodes(i).poly, nodes(j).poly);

            % If pair of nodes overlap, calculate the weight for overlap and record to A
            if ~isempty(ov) && area(ov) > P.overlapThreshold
                w = norm(ci - cj);
                A(i,j)=w; 
                A(j,i)=w;
                continue;
            end
        end
    end
end


function id = findNodeContainingPoint(nodes, q)
    id = [];
    for i=1:numel(nodes)
        if isinterior(nodes(i).poly, q(1), q(2))
            id = i;
            return;
        end
    end
end

function [path, distVal] = dijkstraSparse(A, s, t)
    % Simple Dijkstra for sparse adjacency matrix A (nonnegative weights).
    n = size(A,1);
    dist = inf(n,1); dist(s)=0;
    prev = zeros(n,1);
    visited = false(n,1);

    for iter=1:n
        % pick unvisited node with smallest dist
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

    % reconstruct
    path = t;
    while path(1) ~= s
        path = [prev(path(1)); path]; %#ok<AGROW>
        if path(1)==0
            path = []; distVal=inf; return;
        end
    end
    distVal = dist(t);
end
