function [xdot, U, M, B_prop] = tugboat3d(x, u)
%TUGBOAT3D  3-DOF tugboat model obtained from
% İzzet Kağan Erünsal Thesis
%
% State:
%   x = [eta; nu] = [X; Y; psi; ub; vb; r]
%
% Input:
%   u = [F_L; F_R]
%
% Model:
%   eta_dot = J(psi) * nu
%   M*nu_dot + C(nu)*nu + D(nu)*nu = tau
%   tau = B_prop * input

    % Parameters 
    m  = 10.2;        % [kg] approximate boat mass
    Iz = 0.63994;     % [kg m^2] identified yaw inertia

    % Geometry
    beam = 0.29;      % [m]
    d = beam/2;       % assumed |y_thruster - y_CG| [m]

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
    psi = x(3);
    ub  = x(4);
    vb  = x(5);
    r   = x(6);

    FL = u(1);
    FR = u(2);

    eta = x(1:3);
    nu  = x(4:6); %nu = [ub; vb; r];
    U  = hypot(ub, vb);

    % Thruster allocation: tau = [X_force; Y_force; N_moment]

    B_prop = [1, 1;
              0, 0;
              d, -d];

    tau = B_prop * [FL; FR];

    % 3-DOF inertia matrix
    M = [m - Xudot,      0,            0;
         0,              m - Yvdot,   -Yrdot;
         0,             -Yrdot,        Iz - Nrdot];

    % Coriolis and centripetal matrix
    m11 = M(1,1);
    m22 = M(2,2);
    m23 = M(2,3);

    C = [0,                0,  -(m22*vb + m23*r);
         0,                0,              m11*ub;
         m22*vb + m23*r,  -m11*ub,              0];

    % Linear and nonlinear damping matrix
    D = [-Xu - Xuu*abs(ub),                  0,                 0;
                          0, -Yv - Yvv*abs(vb),               -Yr;
                          0,                -Nv, -Nr - Nrr*abs(r)];



    % Dynamics
    nu_dot = M \ (tau - C*nu - D*nu);

    % Kinematics
    J = [cos(psi), -sin(psi), 0;
         sin(psi),  cos(psi), 0;
            0,         0, 1];

    eta_dot = J * nu;

    % Complete state derivative
    xdot = [eta_dot;
            nu_dot];
end


