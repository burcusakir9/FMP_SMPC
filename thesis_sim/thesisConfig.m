function cfg = thesisConfig(vesselName, modelOpts)
%THESISCONFIG  Common simulation setup for all thesis experiments.
%
%   cfg = thesisConfig()                         tugboat (default)
%   cfg = thesisConfig('otter')                  Otter USV
%   cfg = thesisConfig('tugboat', modelOpts)     tugboat with overridden assumed
%                                                values (modelOpts.Fmin, .T_act)
%
% Every experiment uses these values; an experiment only changes the speed
% law, the mission speed, the safety filter, the disturbance and, for the
% sensitivity study, the vessel assumptions and T_R.

    if nargin < 1, vesselName = 'tugboat'; end
    if nargin < 2, modelOpts = struct(); end

    here = fileparts(mfilename('fullpath'));
    addpath(fileparts(here));                   % vessel models, getScenario, map.kml

    %% Map and RSC funnel chains
    cfg.scenario = 6;                           % tug-scale harbour, 120 x 80 m (getScenario.m)
    cfg.seeds    = 1:2;                         % one RSC chain per seed (supporting
                                                % simulations; use e.g. 1:20 for statistics)
    if ~isempty(getenv('THESIS_SEEDS'))         % quick test: setenv('THESIS_SEEDS', '1:2')
        cfg.seeds = eval(getenv('THESIS_SEEDS'));
    end
    cfg.rsc = struct( ...
        'minCircArea',  3.0, ...                % [m^2] minimum funnel area
        'maxRadius',    inf, ...                % [m]
        'safetyMargin', 2.0, ...                % [m] clearance to obstacles
        'expandStep',   0.1, ...                % [m]
        'eta',          0.9, ...                % child center at eta*R_parent
        'maxIter',      500);

    %% Vessel
    switch vesselName
        case 'tugboat'                          % Erunsal (2015) tug model, 26 N per thruster
            cfg.model  = @(x, u) tugboat3d(x, u, modelOpts);
            cfg.vessel = tugboat3d([], [], modelOpts);
        case 'otter'                            % Fossen's Otter USV, 120 N per propeller
            cfg.model  = @otter3d;
            cfg.vessel = otter3d();
        otherwise
            error('thesisConfig: unknown vessel "%s"', vesselName);
    end
    cfg.vesselName = vesselName;

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
    %   'durmaz' : s = 2*Kv*rho in every funnel (original, unchanged)
    %   'ms'     : s = U_m in intermediate funnels (mission speed)
    %   'msR'    : s = U_k = min(U_m, R_k/T_R) in intermediate funnels
    %   goal funnel of 'ms' and 'msR': s = max(U*tanh(2*Kv*rho/U), min(U, U_goalMin))
    %   (U = U_m or U_k), so the vessel never slows below U_goalMin before it
    %   reaches the goal tolerance and can still close against a current
    cfg.T_R       = 1/cfg.Ka;                   % [s] heading time constant
    cfg.U_goalMin = 0.5;                        % [m/s] minimum speed in the goal funnel

    %% Low-level PI loops (lowLevelControl.m tugboat gains, scaled with the
    %% vessel's surge mass and yaw inertia so the loop bandwidths are equal)
    s_u = cfg.vessel.M(1,1)/11.07818;
    s_r = cfg.vessel.M(3,3)/0.68525;
    cfg.PI = struct('P_speed', 100*s_u, 'I_speed', 50*s_u, 'P_yaw', 5*s_r, 'I_yaw', 0.02*s_r, ...
                    'Ispeed_max', 2000*cfg.vessel.Fmax, 'Iyaw_max', 20000);

    %% Safety filters (QPs with a heavily penalized slack, so they always
    %% return the input that violates the barrier condition least)
    %   'cbf' / 'hocbf'   nominal barrier conditions
    %   'rcbf' / 'rhocbf' robust: b_dot is replaced by b_dot - Vb, i.e. the
    %                     conditions hold for any current up to Vb [m/s]
    cfg.KF = struct('k1', 5, 'k2', 5, 'w_lim', pi/2, 'slack_w', 1e4);   % kinematic CBF/HOCBF on (u, w)
    cfg.HF = struct('k1', 5, 'k2', 5, 'slack_w', 1e4);                  % dynamic HOCBF on thrust
    cfg.Vb = 0.3;                               % [m/s] current bound of the robust filters
    cfg.qp = optimoptions('quadprog', 'Display', 'off');

    %% Output folders
    cfg.dirCache   = fullfile(here, 'cache');
    cfg.dirResults = fullfile(here, 'results');
    cfg.dirFigures = fullfile(here, 'figures');
    for d = {cfg.dirCache, cfg.dirResults, cfg.dirFigures}
        if ~exist(d{1}, 'dir'), mkdir(d{1}); end
    end
end
