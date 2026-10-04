function xdot = tugboat3d(x, u, opts)
%TUGBOAT3D  3-DOF tugboat model with first-order actuator dynamics
%
%   xdot = tugboat3d(x, u)   state derivative
%   p    = tugboat3d()       parameters and actuator limits (struct)
%   ...  = tugboat3d(x, u, opts) / tugboat3d([], [], opts)
%          opts.Fmin, opts.T_act override the two assumed values (sensitivity)
%
% Vessel: 1/40 Pacific Islander tug model of Erunsal (2015), "System
% identification and control of a sea surface vehicle", MSc thesis, METU.
% 900 mm long, 290 mm beam, 10.2 kg without batteries (Sec. 5.2.2), two
% independent stern propellers in Kort nozzles. Hydrodynamic derivatives:
% thesis Table 3.8.
%
% State:
%   x = [eta; nu; F]
%     = [X; Y; psi; ub; vb; r; F_L; F_R]
%
% Input:
%   u = [F_L_cmd; F_R_cmd]   saturated to [Fmin, Fmax] per thruster
%
% Model:
%   eta_dot = J(psi)*nu
%   M*nu_dot + C(nu)*nu + D(nu)*nu = tau
%   tau = B_prop*F
%
% Actuator model:
%   T_L*F_L_dot + F_L = F_L_cmd
%   T_R*F_R_dot + F_R = F_R_cmd
%
% Limits:
%   Fmax = 26 N   largest thrust measured per thruster (thesis Fig. 3.13)
%   Fmin = -14.5 N  ASSUMED: reverse thrust was not measured; set to the
%                   Otter's reverse/forward bollard ratio (13.6/24.4)
%   Umax = 2 m/s  reported top speed (Erunsal et al., ICCAS 2017)
%   T_L, T_R      ASSUMED: thruster dynamics were not identified

    % Parameters
    m  = 10.2;        % [kg] boat mass (thesis Sec. 5.2.2)
    Iz = 0.63994;     % [kg m^2] identified yaw inertia

    % Geometry
    L    = 0.90;      % [m]
    beam = 0.29;      % [m]
    d = beam/2;       % assumed |y_thruster - y_CG| [m]

    % Actuator limits and time constants
    Fmax = 26.0;                  % [N] forward thrust per thruster
    Fmin = -Fmax*13.6/24.4;       % [N] reverse thrust per thruster (assumed)
    Umax = 2.0;                   % [m/s] reported top speed
    T_L  = 0.25;                  % [s] left actuator time constant (assumed)
    T_R  = 0.25;                  % [s] right actuator time constant (assumed)
    if nargin < 3, opts = struct(); end
    if isfield(opts, 'Fmin'),  Fmin = opts.Fmin; end
    if isfield(opts, 'T_act'), T_L = opts.T_act; T_R = opts.T_act; end

    % Identified hydrodynamic derivatives
    Xu    = -5.76909;
    Xuu   = -2.17161;
    Xudot = -0.87818;

    Yv    = -3.98659;
    Yr    = -0.0001;
    Yvv   = -3.95131;
    Yvdot = -1.05279;
    Yrdot = -1.92760;

    Nv    = -0.0001;
    Nr    = -0.12392;
    Nrr   = -0.33077;
    Nrdot = -0.04531;

    % 3-DOF inertia matrix
    M = [m - Xudot,      0,            0;
         0,              m - Yvdot,   -Yrdot;
         0,             -Yrdot,        Iz - Nrdot];

    % Thruster allocation
    B_prop = [1, 1;
              0, 0;
              d, -d];

    if nargin == 0 || isempty(x)
        xdot = struct('name', 'tugboat', 'L', L, 'beam', beam, 'm', m, 'd', d, ...
            'Fmax', Fmax, 'Fmin', Fmin, 'Umax', Umax, 'T_act', T_L, 'M', M, 'B_prop', B_prop);
        return;
    end

    % States
    eta = x(1:3);
    nu  = x(4:6);
    F   = x(7:8);

    psi = eta(3);

    ub = nu(1);
    vb = nu(2);
    r  = nu(3);

    FL = F(1);
    FR = F(2);

    % Commanded actuator forces, saturated to the thruster limits
    FL_cmd = min(max(u(1), Fmin), Fmax);
    FR_cmd = min(max(u(2), Fmin), Fmax);

    % Actual actuator forces are applied to the vessel
    tau = B_prop * F;

    % Coriolis and centripetal matrix
    m11 = M(1,1);
    m22 = M(2,2);
    m23 = M(2,3);

    C = [0,                0, -(m22*vb + m23*r);
         0,                0,             m11*ub;
         m22*vb + m23*r, -m11*ub,               0];

    % Linear and nonlinear damping matrix
    D = [-Xu - Xuu*abs(ub),                  0,                 0;
                          0, -Yv - Yvv*abs(vb),               -Yr;
                          0,                -Nv, -Nr - Nrr*abs(r)];

    % Vessel dynamics
    nu_dot = M \ (tau - C*nu - D*nu);

    % Kinematics
    J = [cos(psi), -sin(psi), 0;
         sin(psi),  cos(psi), 0;
                0,         0, 1];

    eta_dot = J * nu;

    % First-order actuator dynamics
    FL_dot = (FL_cmd - FL) / T_L;
    FR_dot = (FR_cmd - FR) / T_R;

    F_dot = [FL_dot;
             FR_dot];

    % Complete state derivative
    xdot = [eta_dot;
            nu_dot;
            F_dot];
end
