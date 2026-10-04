function cfg = thesisConfig()
%THESISCONFIG  Common simulation setup for all thesis experiments (H1, H2, H3).
%
% Every experiment uses these values; an experiment only changes the speed
% law, the mission speed, the safety filter and the disturbance of a run.

    here = fileparts(mfilename('fullpath'));
    addpath(fileparts(here));                   % vessel models, getScenario, map.kml

    %% Map and RSC funnel chains
    cfg.scenario = 4;                           % real coastline (getScenario.m, map.kml)
    cfg.seeds    = 1:5;                         % one RSC chain per seed
    cfg.rsc = struct( ...
        'minCircArea',  3.0, ...                % [m^2] minimum funnel area
        'maxRadius',    inf, ...                % [m]
        'safetyMargin', 2.0, ...                % [m] clearance to obstacles
        'expandStep',   0.1, ...                % [m]
        'eta',          0.9, ...                % child center at eta*R_parent
        'maxIter',      500);

    %% Vessel
    cfg.model  = @tugboat3d;                    % Erunsal (2015) tug model, 26 N per thruster
    cfg.vessel = tugboat3d();                   % parameters and thrust limits

    %% Simulation
    cfg.dt       = 0.01;                        % [s] RK4 step
    cfg.goal_tol = 2.0;                         % [m] goal reached
    cfg.Tfactor  = 4;                           % Tmax = Tfactor * path length / U_m (>= 500 s)
    cfg.Tstuck   = 300;                         % [s] stop if no new funnel is entered for this long
    cfg.logEvery = 10;                          % log every 10th step (0.1 s)

    %% Funnel controller (Durmaz2024, Eq. 33)
    cfg.Kv      = 0.05;                         % original law: s = 2*Kv*rho
    cfg.Ka      = 0.30;                         % heading gain: alpha_dot = -Ka*alpha
    cfg.rho_tol = 0.05;                         % [m] stop at the center

    %% Speed laws (u = s(rho)*cos(alpha))
    %   'durmaz' : s = 2*Kv*rho in every funnel (original)
    %   'ms'     : s = U_m in intermediate funnels (mission speed)
    %   'msR'    : s = U_k = min(U_m, R_k/T_R) in intermediate funnels
    %              (mission speed scheduled by funnel radius)
    %   goal funnel of 'ms' and 'msR': s = U*tanh(2*Kv*rho/U), U = U_m or U_k
    cfg.T_R = 1/cfg.Ka;                         % [s] heading time constant

    %% Low-level PI loops (lowLevelControl.m gains) + allocation
    cfg.PI = struct('P_speed', 100, 'I_speed', 50, 'P_yaw', 5, 'I_yaw', 0.02, ...
                    'Ispeed_max', 2000*cfg.vessel.Fmax, 'Iyaw_max', 20000);

    %% Safety filters (both QPs with a heavily penalized slack, so they always
    %% return the input that violates the barrier condition least)
    cfg.KF = struct('k1', 5, 'k2', 5, 'w_lim', pi/2, 'slack_w', 1e4);   % kinematic CBF/HOCBF on (u, w)
    cfg.HF = struct('k1', 5, 'k2', 5, 'slack_w', 1e4);                  % dynamic HOCBF on thrust
    cfg.qp = optimoptions('quadprog', 'Display', 'off');

    %% Output folders
    cfg.dirCache   = fullfile(here, 'cache');
    cfg.dirResults = fullfile(here, 'results');
    cfg.dirFigures = fullfile(here, 'figures');
    for d = {cfg.dirCache, cfg.dirResults, cfg.dirFigures}
        if ~exist(d{1}, 'dir'), mkdir(d{1}); end
    end
end
