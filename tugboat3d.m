function xdot = tugboat3d(x, u)
%TUGBOAT3D  3-DOF tugboat model with first-order actuator dynamics
%
% State:
%   x = [eta; nu; F]
%     = [X; Y; psi; ub; vb; r; F_L; F_R]
%
% Input:
%   u = [F_L_cmd; F_R_cmd]
%
% Model:
%   eta_dot = J(psi)*nu
%   M*nu_dot + C(nu)*nu + D(nu)*nu = tau
%   tau = B_prop*F
%
% Actuator model:
%   T_L*F_L_dot + F_L = F_L_cmd
%   T_R*F_R_dot + F_R = F_R_cmd

    % Parameters
    m  = 10.2;        % [kg] approximate boat mass
    Iz = 0.63994;     % [kg m^2] identified yaw inertia

    % Geometry
    beam = 0.29;      % [m]
    d = beam/2;       % assumed |y_thruster - y_CG| [m]

    % Actuator time constants
    T_L = 0.25;       % [s] left actuator time constant
    T_R = 0.25;       % [s] right actuator time constant

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

    % Commanded actuator forces
    FL_cmd = u(1);
    FR_cmd = u(2);

    U = hypot(ub, vb);

    % Thruster allocation
    B_prop = [1, 1;
              0, 0;
              d, -d];

    % Actual actuator forces are applied to the vessel
    tau = B_prop * F;

    % 3-DOF inertia matrix
    M = [m - Xudot,      0,            0;
         0,              m - Yvdot,   -Yrdot;
         0,             -Yrdot,        Iz - Nrdot];

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