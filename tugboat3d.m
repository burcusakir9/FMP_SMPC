function [xdot, U, M, B_prop] = tugboat3d(x, u)
%TUGBOAT3D  3-DOF tugboat model reconstructed from Erunsal/Kumru theses
%
% State:
%   x = [X; Y; psi; u; v; r]
%       X, Y  : inertial position [m]
%       psi   : yaw angle [rad]
%       u     : surge speed in body frame [m/s]
%       v     : sway speed in body frame [m/s]
%       r     : yaw rate [rad/s]
%
% Input:
%   u = [F_L; F_R]
%       F_L   : left thruster force [N]
%       F_R   : right thruster force [N]
%
% Outputs:
%   xdot   : state derivative
%   U      : total planar speed sqrt(u^2 + v^2)
%   M      : 3x3 inertia matrix used in the dynamics
%   B_prop : 3x2 input mapping from [F_L;F_R] to tau=[X;Y;N]
%
% Notes:
% - This is a control-oriented 3-DOF reconstruction. The thesis gives the
%   reduced model symbolically but does not print the full expanded fa...fi terms.
% - Identified parameters are taken from the spiral-maneuver table.
% - Thruster lateral offsets are assumed symmetric at +/- beam/2.
%
% State ordering follows the reduced linearization form [x y psi u v r]^T
% discussed in the thesis.

    %%%%%%%%%%%%%%%
    % Parameters  %
    %%%%%%%%%%%%%%%

    % Basic boat data
    m  = 10.2;        % [kg] approximate boat mass
    Iz = 0.63994;     % [kg m^2] identified yaw inertia

    % Geometry
    beam = 0.29;      % [m]
    d = beam/2;       % assumed |y_thruster - y_CG| [m]

    % 3-DOF identified parameters from spiral-motion results
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

    %%%%%%%%%%%%%%%
    % Unpack data %
    %%%%%%%%%%%%%%%
    psi = x(3);
    ub  = x(4);
    vb  = x(5);
    r   = x(6);

    FL = u(1);
    FR = u(2);

    nu = [ub; vb; r];
    U  = hypot(ub, vb);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Thruster mapping: tau = [X_force; Y_force; N_moment]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % From the thesis, left/right thrusts are mapped into surge force and
    % yaw torque. For symmetric offsets:
    %   tau1 = FL + FR
    %   tau3 = d*(FL - FR)
    %
    B_prop = [1, 1;
              0, 0;
              d, -d];

    tau = B_prop * [FL; FR];

    %%%%%%%%%%%%%%%%%%%%%%%
    % 3-DOF inertia matrix
    %%%%%%%%%%%%%%%%%%%%%%%
    %
    % Added-mass coupling in sway-yaw is retained using Yrdot.
    % This is a compact Fossen-style 3-DOF approximation.
    %
    M = [m - Xudot,      0,            0;
         0,              m - Yvdot,   -Yrdot;
         0,             -Yrdot,        Iz - Nrdot];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Coriolis/centripetal matrix (control-oriented form)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m11 = M(1,1);
    m22 = M(2,2);
    m23 = M(2,3);

    C = [0,  0, -(m22*vb + m23*r);
         0,  0,  m11*ub;
         m22*vb + m23*r, -m11*ub, 0];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Hydrodynamic damping / restoring in 3DOF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Coefficients are already identified with sign, so we apply them
    % directly as generalized forces/moments.
    %
    tau_h = [ Xu*ub + Xuu*abs(ub)*ub;
              Yv*vb + Yr*r + Yvv*abs(vb)*vb;
              Nv*vb + Nr*r + Nrr*abs(r)*r ];

    %%%%%%%%%%%%%%%%%%%%%%%
    % Solve dynamic model %
    %%%%%%%%%%%%%%%%%%%%%%%
    nudot = M \ (tau + tau_h - C*nu);

    %%%%%%%%%%%%%%
    % Kinematics %
    %%%%%%%%%%%%%%
    Rpsi = [cos(psi), -sin(psi);
            sin(psi),  cos(psi)];

    etadot_xy = Rpsi * [ub; vb];
    psidot = r;

    xdot = [etadot_xy(1);
            etadot_xy(2);
            psidot;
            nudot(1);
            nudot(2);
            nudot(3)];
end