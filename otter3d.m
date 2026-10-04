function xdot = otter3d(x, u)
%OTTER3D  3-DOF Otter USV model with first-order actuator dynamics
%
%   xdot = otter3d(x, u)   state derivative
%   p    = otter3d()       parameters and actuator limits (struct)
%
% Vessel: Maritime Robotics Otter USV, Fossen's otter.m (MSS toolbox),
% reduced to surge-sway-yaw as in otter_mpc_sims/otter/otter3d.m, with
% thrust forces instead of propeller speeds as inputs. 2.0 m long, 1.08 m
% beam, 55 kg (no payload), two propellers on the pontoons.
% Same interface as tugboat3d.m.
%
% State:
%   x = [X; Y; psi; ub; vb; r; F_L; F_R]
%
% Input:
%   u = [F_L_cmd; F_R_cmd]   saturated to [Fmin, Fmax] per propeller
%
% Model:
%   eta_dot = J(psi)*nu
%   M*nu_dot + C(nu)*nu = tau + tau_damp(nu)
%   tau = B_prop*F,  tau_damp = [Xu*u; Yv*v; Nr*(1 + 10|r|)*r]
%
% Actuator model:
%   T_n*F_dot + F = F_cmd
%
% Limits (Fossen's otter.m bollard pull data):
%   Fmax =  0.5*24.4*g = 119.7 N per propeller (forward)
%   Fmin = -0.5*13.6*g = -66.7 N per propeller (reverse)
%   Umax = 6 knots = 3.09 m/s
%   T_n  = 0.1 s propeller time constant

    % Main data
    g = 9.81;
    L = 2.0;                    % length [m]
    B = 1.08;                   % beam [m]
    m = 55.0;                   % mass [kg]
    rg = [0.2 0 -0.2]';         % CG for hull only [m]
    R44 = 0.4*B;                % radii of gyration [m]
    R55 = 0.25*L;
    R66 = 0.25*L;
    T_sway = 1;                 % time constant in sway [s]
    T_yaw  = 1;                 % time constant in yaw [s]
    Umax = 6*0.5144;            % maximum forward speed [m/s]
    y_pont = 0.395;             % propeller lever arm from centerline [m]

    % Actuator limits and time constant
    Fmax =  0.5*24.4*g;         % [N] forward thrust per propeller
    Fmin = -0.5*13.6*g;         % [N] reverse thrust per propeller
    T_n  = 0.1;                 % [s] propeller time constant

    % Inertia
    Ig = m*diag([R44^2, R55^2, R66^2]) - m*Smtrx(rg)^2;
    Iz = Ig(3,3);

    % Added mass (best practice, Fossen)
    Xudot = -0.1*m;
    Yvdot = -1.5*m;
    Nrdot = -1.7*Iz;

    MRB = diag([m, m, Iz]);
    MA  = -diag([Xudot, Yvdot, Nrdot]);
    M   = MRB + MA;

    % Linear damping
    Xu = -24.4*g/Umax;          % from the maximum speed
    Yv = -M(2,2)/T_sway;        % from the time constant in sway
    Nr = -M(3,3)/T_yaw;         % from the time constant in yaw

    % Thruster allocation (left propeller at +y_pont gives positive yaw)
    B_prop = [1, 1;
              0, 0;
              y_pont, -y_pont];

    if nargin == 0
        xdot = struct('name', 'otter', 'L', L, 'beam', B, 'm', m, 'd', y_pont, ...
            'Fmax', Fmax, 'Fmin', Fmin, 'Umax', Umax, 'T_act', T_n, 'M', M, 'B_prop', B_prop);
        return;
    end

    % States
    psi = x(3);
    nu  = x(4:6);
    F   = x(7:8);

    % Commanded propeller forces, saturated to the bollard pull limits
    F_cmd = min(max(u(:), Fmin), Fmax);

    % Coriolis and centripetal matrices
    CRB = [0,          0,         -m*nu(2);
           0,          0,          m*nu(1);
           m*nu(2),   -m*nu(1),    0];
    CA  = [0,              0,              Yvdot*nu(2);
           0,              0,             -Xudot*nu(1);
          -Yvdot*nu(2),    Xudot*nu(1),    0];
    C = CRB + CA;

    % Linear damping + nonlinear yaw damping
    tau_damp = [Xu*nu(1);
                Yv*nu(2);
                Nr*(1 + 10*abs(nu(3)))*nu(3)];

    % Actual propeller forces are applied to the vessel
    tau = B_prop*F;

    % Kinematics
    J = [cos(psi), -sin(psi), 0;
         sin(psi),  cos(psi), 0;
                0,         0, 1];

    xdot = [J*nu;
            M \ (tau + tau_damp - C*nu);
            (F_cmd - F)/T_n];
end

function S = Smtrx(a)
    S = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
end
