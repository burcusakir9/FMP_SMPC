% [map edges, obstacle polygons, start point, goal point]
function [W, obs, q_start, q_goal] = getScenario(scn)

switch scn
    case 1

        W = [0 30 0 20];

        obs = {};
        obs{end+1} = polyshape([6 10 10 6],[3 3 7 7]); % [x1 x2 x3 x4],[y1 y2 y3 y4]
        obs{end+1} = polyshape([14 18 18 14],[2 2 6 6]);
        obs{end+1} = polyshape([20 24 24 20],[10 10 16 16]);
        obs{end+1} = polyshape([8 12 12 8],[12 12 18 18]);
        obs{end+1} = polyshape([2 4 4 2],[9 9 14 14]);

        q_start = [2 2]; 
        q_goal  = [28 18];

    case 2

        W = [0 16 0 8];
        obs = {};

        obs{end+1} = polyshape([0 1.5 1.5 0], [3.5 3.5 6.0 6.0]);
        obs{end+1} = polyshape([0 6.0 2.3 0], [0 0 2.0 2.0]);
        obs{end+1} = polyshape([3.2 3.2 4.7 7.6 6.9], [3.5 8.0 8.0 5.0 1.2]);
        obs{end+1} = polyshape([7.2 9.6 9.6 7.9], [1.0 0.0 4.5 5.0]);
        obs{end+1} = polyshape([7.5 9.5 9.5 7.5], [6.0 6.5 8.0 8.0]);
        c = [12.8 4.5];   
        r = 1.5;         
        ang = deg2rad(22.5 + (0:7)*45);
        obs{end+1} = polyshape(c(1)+r*cos(ang), c(2)+r*sin(ang));
        obs{end+1} = polyshape([11.0 14.0 14.0 11.0], [0.0 0.0 1.5 1.5]);


        q_start = [2.5 5]; 
        q_goal  = [15 5];
        
    case 3

        W = [0 500 0 500];
        obs = {};
        obs{end+1} = polyshape([200 300 300 200], [50 50 75 75]);
        obs{end+1} = polyshape([225 275 275 225], [225 225 275 275]);
        obs{end+1} = polyshape([375 400 400 375], [150 150 300 300]);
        obs{end+1} = polyshape([50 75 75 50], [275 275 325 325]);
        obs{end+1} = polyshape([400 450 450 400], [25 25 75 75]);
        q_start = [25 25]; 
        q_goal  = [475 450];


    case 4 % Real KML Coastline: Cropped, Scaled (300x300), and Padded
        try
            kmlData = readstruct('map.kml', 'FileType', 'xml');
        catch
            error('map.kml not found. Ensure it is in the current directory.');
        end

        placemarks = kmlData.Document.Placemark;
        numPolys = length(placemarks);
        raw_obs = {};
        all_x_raw = []; all_y_raw = [];

        % Reference Origin (First point of first polygon)
        firstCoordStr = strtrim(placemarks(1).Polygon.outerBoundaryIs.LinearRing.coordinates);
        firstPointSplit = split(firstCoordStr, ' ');
        firstPointCoords = split(firstPointSplit{1}, ',');
        lon0 = str2double(firstPointCoords{1});
        lat0 = str2double(firstPointCoords{2});

        % Crop & Conversion Constants
        y_max_limit = 3300; 
        clipper = polyshape([-1e7 1e7 1e7 -1e7], [-1e7 -1e7 y_max_limit y_max_limit]);

        % First Pass: Meter Conversion and Cropping
        for i = 1:numPolys
            coordStr = regexprep(strtrim(placemarks(i).Polygon.outerBoundaryIs.LinearRing.coordinates), '\s+', ' ');
            points = split(coordStr, ' ');
            curLat = []; curLon = [];
            for j = 1:length(points)
                if isempty(points{j}), continue; end
                val = split(points{j}, ',');
                if length(val) >= 2
                    curLon(end+1) = str2double(val{1});
                    curLat(end+1) = str2double(val{2});
                end
            end
            x = 111320 * (curLon - lon0) * cosd(lat0);
            y = 110574 * (curLat - lat0);
            pgon_c = intersect(polyshape(x, y), clipper);
            if pgon_c.NumRegions > 0
                raw_obs{end+1} = pgon_c;
                [xc, yc] = boundary(pgon_c);
                all_x_raw = [all_x_raw; xc]; all_y_raw = [all_y_raw; yc];
            end
        end

        % Second Pass: Scaling (300x300) and Padding (L/R: 50, B: 100)
        minX = min(all_x_raw); maxX = max(all_x_raw);
        minY = min(all_y_raw); maxY = max(all_y_raw);
        scaleX = 300 / (maxX - minX);
        scaleY = 300 / (maxY - minY);

        obs = {};
        for i = 1:length(raw_obs)
            [vx, vy] = boundary(raw_obs{i});
            new_vx = ((vx - minX) * scaleX) + 50;  % 50 padding left
            new_vy = ((vy - minY) * scaleY) + 100; % 100 padding bottom
            obs{end+1} = polyshape(new_vx, new_vy);
        end

        % --- ADD DETAILED IRREGULAR ISLANDS (10+ Edges Each) ---
        
        % Island 1: Jagged "C" Rock (Bottom-Left Buffer)
        % 12 vertices: Creates a complex shoreline with a small bay
        v1_x = [70 85 100 115 130 125 110 95 80 65 60 55];
        v1_y = [15 10 15 30 50 65 60 45 60 50 35 25];
        obs{end+1} = polyshape(v1_x, v1_y);
        
        
        % Island 3: Complex Central Hazard (Center-Left)
        % 10 vertices: A bulky, uneven mass
        v3_x = [35 65 85 80 60 45 25 15 20 30];
        v3_y = [220 215 235 255 270 265 255 240 225 215];
        obs{end+1} = polyshape(v3_x, v3_y);
        
        W = [0 400 0 400]; % Total space: 50+300+50 by 100+300
        q_start = [50 350];  % Center-bottom in the 100-unit padding
        q_goal  = [350 350]; % Inside the upper mapped area


    otherwise
        error('Unknown scenario id: %d', scn);
end
end
