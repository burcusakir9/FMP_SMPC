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

    otherwise
        error('Unknown scenario id: %d', scn);
end
end
