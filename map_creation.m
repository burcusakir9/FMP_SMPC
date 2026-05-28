clear; clc; close all;

% 1. Load the KML file
try
    kmlData = readstruct('map.kml', 'FileType', 'xml');
catch
    error('File "map.kml" not found. Ensure it is in the current directory.');
end

% 2. Extract Placemarks
placemarks = kmlData.Document.Placemark;
numPolys = length(placemarks);
raw_obs = {};
all_x_raw = [];
all_y_raw = [];

% 3. Set Reference Origin
firstCoordStr = strtrim(placemarks(1).Polygon.outerBoundaryIs.LinearRing.coordinates);
firstPointSplit = split(firstCoordStr, ' ');
firstPointCoords = split(firstPointSplit{1}, ',');
lon0 = str2double(firstPointCoords{1});
lat0 = str2double(firstPointCoords{2});

% 4. Crop Geometry (at y = 3300 meters)
y_max_limit = 3300; 
clipper = polyshape([-1e7 1e7 1e7 -1e7], [-1e7 -1e7 y_max_limit y_max_limit]);

% 5. First Pass: Convert and Crop
for i = 1:numPolys
    coordStr = strtrim(placemarks(i).Polygon.outerBoundaryIs.LinearRing.coordinates);
    coordStr = regexprep(coordStr, '\s+', ' '); 
    points = split(coordStr, ' ');
    
    currentLat = []; currentLon = [];
    for j = 1:length(points)
        if isempty(points{j}), continue; end
        val = split(points{j}, ',');
        if length(val) >= 2
            currentLon(end+1) = str2double(val{1});
            currentLat(end+1) = str2double(val{2});
        end
    end
    
    x = 111320 * (currentLon - lon0) * cosd(lat0);
    y = 110574 * (currentLat - lat0);
    
    pgon_cropped = intersect(polyshape(x, y), clipper);
    
    if pgon_cropped.NumRegions > 0
        raw_obs{end+1} = pgon_cropped;
        [xc, yc] = boundary(pgon_cropped);
        all_x_raw = [all_x_raw; xc];
        all_y_raw = [all_y_raw; yc];
    end
end

% 6. Second Pass: Scale to 300x300 and Add Padding
% Padding Definitions
pad_left = 50;
pad_right = 50;
pad_bottom = 100;
target_size = 300;

% Current raw bounds
minX = min(all_x_raw); maxX = max(all_x_raw);
minY = min(all_y_raw); maxY = max(all_y_raw);

% Scale factors to fit the 300x300 box
scaleX = target_size / (maxX - minX);
scaleY = target_size / (maxY - minY);

obs = {};
for i = 1:length(raw_obs)
    [vx, vy] = boundary(raw_obs{i});
    
    % Scale to 300x300, then shift by padding
    % X: (Normalize * Scale) + 50 (left pad)
    % Y: (Normalize * Scale) + 100 (bottom pad)
    new_vx = ((vx - minX) * scaleX) + pad_left;
    new_vy = ((vy - minY) * scaleY) + pad_bottom;
    
    obs{end+1} = polyshape(new_vx, new_vy);
end

% 7. Define Final Workspace
% X-axis: 50 (left) + 300 (map) + 50 (right) = 400
% Y-axis: 100 (bottom) + 300 (map) = 400
W = [0 400 0 400];
q_start = [20, 20]; % Inside the 100-unit bottom buffer
q_goal  = [200, 350]; % Inside the mapped area

% --- Visualization ---
figure('Color', 'w', 'Name', 'Final Padded Scenario');
hold on;

% Draw the 300x300 "Active Map" boundary for reference
rectangle('Position', [50, 100, 300, 300], 'EdgeColor', [0.8 0.8 0.8], 'LineStyle', '--');

% Plot Padded Polygons
for i = 1:length(obs)
    plot(obs{i}, 'FaceColor', [0.2 0.6 1.0], 'FaceAlpha', 0.5, 'EdgeColor', 'k');
end

% Plot Start/Goal
plot(q_start(1), q_start(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
plot(q_goal(1), q_goal(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);

% Final Formatting
axis equal; grid on; box on;
xlim([0 400]); ylim([0 400]);
xlabel('X (meters)'); ylabel('Y (meters)');
title('KML: 300x300 Map with Padding (L:50, R:50, B:100)');
legend('Map Boundary', 'Obstacles', 'Start', 'Goal', 'Location', 'northeastoutside');

fprintf('Successfully scaled map to 300x300 and applied padding.\n');
fprintf('Final Workspace: [0, 400, 0, 400]\n');