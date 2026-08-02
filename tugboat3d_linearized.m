function [A_c, B_c] = tugboat3d_linearized(x_bar, u_bar)
%TUGBOAT3D_LINEARIZED Continuous-time analytical linearization of tugboat3d
%
% Linearization point:
%   x_bar = [eta_bar; nu_bar; F_bar]
%         = [X; Y; psi; ub; vb; r; F_L; F_R]
%
%   u_bar = [F_L_cmd; F_R_cmd]
%
% Linear error-state model:
%   delta_xdot = A_c*delta_x + B_c*delta_u

    % Ensure column vectors
    x_bar = x_bar(:);
    u_bar = u_bar(:);

    % Validate dimensions
    if numel(x_bar) ~= 8
        error('x_bar must contain 8 elements.');
    end

    if numel(u_bar) ~= 2
        error('u_bar must contain 2 elements.');
    end

    % Parameters
    m  = 10.2;
    Iz = 0.63994;

    % Geometry
    beam = 0.29;
    d = beam/2;

    % Actuator time constants
    T_L = 0.25;
    T_R = 0.25;

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

    % Nominal states
    eta_bar = x_bar(1:3);
    nu_bar  = x_bar(4:6);
    F_bar   = x_bar(7:8);

    psi = eta_bar(3);

    ub = nu_bar(1);
    vb = nu_bar(2);
    r  = nu_bar(3);

    FL = F_bar(1);
    FR = F_bar(2);

    % Thruster allocation
    B_prop = [1, 1;
              0, 0;
              d, -d];

    % 3-DOF inertia matrix
    M = [m - Xudot,      0,            0;
         0,              m - Yvdot,   -Yrdot;
         0,             -Yrdot,        Iz - Nrdot];

    % Continuous-time state Jacobian
    A_c = zeros(8,8);

    cpsi = cos(psi);
    spsi = sin(psi);

    % Kinematic Jacobian
    %
    % eta_dot = J(psi)*nu
    A_c(1,3) = -spsi*ub - cpsi*vb;
    A_c(1,4) =  cpsi;
    A_c(1,5) = -spsi;

    A_c(2,3) =  cpsi*ub - spsi*vb;
    A_c(2,4) =  spsi;
    A_c(2,5) =  cpsi;

    A_c(3,6) = 1;

    % Dynamic Jacobian
    %
    % nu_dot = M \ (B_prop*F - C(nu)*nu - D(nu)*nu)
    m11 = M(1,1);
    m22 = M(2,2);
    m23 = M(2,3);

    dF_dnu = zeros(3,3);

    % Surge equation derivatives
    dF_dnu(1,1) = Xu + 2*Xuu*abs(ub);
    dF_dnu(1,2) = m22*r;
    dF_dnu(1,3) = m22*vb + 2*m23*r;

    % Sway equation derivatives
    dF_dnu(2,1) = -m11*r;
    dF_dnu(2,2) = Yv + 2*Yvv*abs(vb);
    dF_dnu(2,3) = Yr - m11*ub;

    % Yaw equation derivatives
    dF_dnu(3,1) = -(m22 - m11)*vb - m23*r;
    dF_dnu(3,2) = Nv - (m22 - m11)*ub;
    dF_dnu(3,3) = Nr + 2*Nrr*abs(r) - m23*ub;

    A_c(4:6,4:6) = M \ dF_dnu;

    % Effect of actual actuator forces on vessel dynamics
    A_c(4:6,7:8) = M \ B_prop;

    % First-order actuator state dynamics
    %
    % F_L_dot = (F_L_cmd - F_L)/T_L
    % F_R_dot = (F_R_cmd - F_R)/T_R
    A_c(7,7) = -1/T_L;
    A_c(8,8) = -1/T_R;

    % Continuous-time input Jacobian
    B_c = zeros(8,2);

    % Commanded forces affect the actuator states
    B_c(7,1) = 1/T_L;
    B_c(8,2) = 1/T_R;
end