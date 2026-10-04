function xdot = cybership3d(x, u)
%CYBERSHIP3D  3-DOF CyberShip II model with first-order actuator dynamics
%
%   xdot = cybership3d(x, u)   state derivative
%   p    = cybership3d()       parameters and actuator limits (struct)
%
% Vessel: CyberShip II, 1:70 supply ship model (NTNU), 23.8 kg, 1.255 m.
% Model and parameters from Skjetne, Smogeli & Fossen (2004), "A nonlinear
% ship manoeuvering model: identification and adaptive control with
% experiments for a model ship", MIC 25(1): Table 1 (mass, added mass,
% actuator positions), Table 2 (surge/sway damping, thrust coefficients),
% Table 4 (yaw-related damping). All taken with respect to CP.
% Same interface as tugboat3d.m.
%
% ACTUATION: only the two main propellers are modelled, as a differential
% thrust pair at (lx, ly) = (-0.499, -/+0.078) m. The real ship also has two
% rudders and a bow thruster, which are left out so that all three vessels
% share the input u = [F_L; F_R]. With a lever arm of only 0.078 m, yaw
% authority is much weaker than on the real ship.
%
% State:
%   x = [X; Y; psi; ub; vb; r; F_L; F_R]
%
% Input:
%   u = [F_L_cmd; F_R_cmd]   saturated to [Fmin, Fmax] per propeller
%
% Model:
%   eta_dot = J(psi)*nu
%   M*nu_dot + C(nu)*nu + D(nu)*nu = tau,   tau = B_prop*F
%
% Limits:
%   Fmax =  T+|n|n * n^2 =  4.1 N   at n = 2000 rpm, the highest propeller
%   Fmin = -T-|n|n * n^2 = -5.7 N   speed in the towing tests (bollard values;
%                                   the paper gives no hardware rpm limit)
%   Umax ~ 1.0 m/s                  steady speed at 2*Fmax (from the damping)
%   T_act = 0.25 s                  ASSUMED (not identified), as tugboat3d.m

    % Mass-related parameters (Table 1)
    m     = 23.8;       % [kg]
    Iz    = 1.760;      % [kg m^2]
    xg    = 0.046;      % [m]
    Xudot = -2.0;
    Yvdot = -10.0;
    Yrdot = -0.0;
    Nvdot = -0.0;
    Nrdot = -1.0;

    % Surge and sway damping (Table 2)
    Xu   = -0.72253;
    Xuu  = -1.32742;    % X_|u|u
    Xuuu = -5.86643;
    Yv   = -0.88965;
    Yvv  = -36.47287;   % Y_|v|v
    Nv   =  0.03130;
    Nvv  =  3.95645;    % N_|v|v

    % Yaw-related damping (Table 4)
    Yrv = -0.805;       % Y_|r|v
    Yr  = -7.250;
    Yvr = -0.845;       % Y_|v|r
    Yrr = -3.450;       % Y_|r|r
    Nrv =  0.130;       % N_|r|v
    Nr  = -1.900;
    Nvr =  0.080;       % N_|v|r
    Nrr = -0.750;       % N_|r|r

    % Geometry and actuators
    L  = 1.255;         % [m]
    B  = 0.29;          % [m]
    lx = -0.499;        % main propellers, x position w.r.t. CP [m]
    ly =  0.078;        % main propellers, |y| position [m]
    n_max = 2000/60;    % [rps] highest tested propeller speed
    Tnn_pos = 3.65034e-3;       % T+_|n|n (Table 2)
    Tnn_neg = 5.10256e-3;       % T-_|n|n (Table 2)
    Fmax =  Tnn_pos*n_max^2;    % [N]
    Fmin = -Tnn_neg*n_max^2;    % [N]
    T_act = 0.25;               % [s] actuator time constant (assumed)

    % Inertia matrix (symmetric, Yrdot = Nvdot = 0)
    M = [m - Xudot,  0,             0;
         0,          m - Yvdot,     m*xg - Yrdot;
         0,          m*xg - Nvdot,  Iz - Nrdot];

    % Thruster allocation: port propeller (y = -ly) gives positive yaw
    B_prop = [1,   1;
              0,   0;
              ly, -ly];

    if nargin == 0
        Umax = fzero(@(s) 2*Fmax + Xu*s + Xuu*abs(s)*s + Xuuu*s^3, [0 5]);
        xdot = struct('name', 'cybership', 'L', L, 'beam', B, 'm', m, 'd', ly, ...
            'Fmax', Fmax, 'Fmin', Fmin, 'Umax', Umax, 'T_act', T_act, 'M', M, 'B_prop', B_prop, ...
            'lx', lx);
        return;
    end

    % States
    psi = x(3);
    nu  = x(4:6);
    F   = x(7:8);
    ub = nu(1); vb = nu(2); r = nu(3);

    % Commanded propeller forces, saturated
    F_cmd = min(max(u(:), Fmin), Fmax);

    % Coriolis and centripetal matrix
    m11 = M(1,1); m22 = M(2,2); m23 = M(2,3);
    C = [0,                0, -(m22*vb + m23*r);
         0,                0,             m11*ub;
         m22*vb + m23*r, -m11*ub,               0];

    % Nonlinear damping matrix
    D = [-Xu - Xuu*abs(ub) - Xuuu*ub^2, 0,                                   0;
         0, -Yv - Yvv*abs(vb) - Yrv*abs(r), -Yr - Yvr*abs(vb) - Yrr*abs(r);
         0, -Nv - Nvv*abs(vb) - Nrv*abs(r), -Nr - Nvr*abs(vb) - Nrr*abs(r)];

    % Actual propeller forces are applied to the vessel
    tau = B_prop*F;

    % Kinematics
    J = [cos(psi), -sin(psi), 0;
         sin(psi),  cos(psi), 0;
                0,         0, 1];

    xdot = [J*nu;
            M \ (tau - C*nu - D*nu);
            (F_cmd - F)/T_act];
end
