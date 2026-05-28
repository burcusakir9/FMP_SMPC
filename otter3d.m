function [xdot, U, M, B_prop] = otter3d(x,n)
% Compatibel with MATLAB and the free software GNU Octave (www.octave.org)
% [xdot,U] = otter(x,n,mp,rp,V_c,beta_c) returns the speed U in m/s 
% (optionally), the 6x6 mass matrix M (optionally), and the 2x4 input 
% matrix B_prop (optionally) in pitch and yaw and the time derivative of 
% the state vector: 
% A
% x = [ x y yaw u v r ]' 
% 
% for the Maritime Robotics Otter USV, see www.maritimerobotics.com. 
% The length of the USV is L = 2.0 m, with the state vector defined as:
%
%  x:     position in x direction (m)
%  y:     position in y direction (m)
%  yaw:   yaw angle               (rad)
%  u:     surge velocity          (m/s)
%  v:     sway velocity           (m/s)
%  r:     yaw velocity            (rad/s)

% Author:    Thor I. Fossen
% Date:      2019-07-17
% Modified:  [Your Name]
% Date:     [Today's Date] - Reordered states to [x, y, yaw, u, v, r]

% Main data
g   = 9.81;         % acceleration of gravity (m/s^2)
rho = 1025;         % density of water
L = 2.0;            % length (m)
B = 1.08;           % beam (m)
m = 55.0;           % mass (kg)
rg = [0.2 0 -0.2]'; % CG for hull only (m)
R44 = 0.4 * B;      % radii of gyrations (m)
R55 = 0.25 * L;
R66 = 0.25 * L;
T_sway = 1;         % time constant in sway (s)
T_yaw = 1;          % time constant in yaw (s)
Umax = 6 * 0.5144;  % 6 knots maximum forward speed (m/s)

% Data for one pontoon
y_pont  = 0.395;    % distance from centerline to waterline area center (m)

% State and current variables
eta = x(1:3);       % position and orientation [x, y, yaw]
nu = x(4:6);        % velocities [u, v, r]
U = sqrt(nu(1)^2 + nu(2)^2 );  % speed

% Inertia dyadic, volume displacement and draft
Ig_CG = m * diag([R44^2, R55^2, R66^2]);    % only hull in the CG
Ig = Ig_CG - m * Smtrx(rg)^2; % hull + payload in the CO

% Experimental propeller data including lever arms
l1 = -y_pont;                           % lever arm, left propeller (m)
l2 = y_pont;                            % lever arm, right propeller (m)
k_pos = 0.02216/2;                      % Positive Bollard, one propeller 
k_neg = 0.01289/2;                      % Negative Bollard, one propeller 
n_max =  sqrt((0.5*24.4 * g)/k_pos);    % maximum propeller rev. (rad/s)
n_min = -sqrt((0.5*13.6 * g)/k_neg);    % minimum propeller rev. (rad/s)

% MRB and CRB (Fossen 2021)
Iz = Ig(3,3);

MRB = [m, 0, 0;
       0, m, 0;
       0, 0, Iz];

CRB = [0      0      -m*nu(2);
       0      0       m*nu(1);
       m*nu(2) -m*nu(1)  0];

% Hydrodynamic added mass (best practice)
Xudot = -0.1 * m;   % -addedMassSurge(m,L,rho);   
Yvdot = -1.5 * m;
Nrdot = -1.7 * Ig(3,3);

MA = -diag([Xudot, Yvdot, Nrdot]);   
CA  = [ 0               0               Yvdot * nu(2);
        0               0              -Xudot * nu(1);
       -Yvdot * nu(2)   Xudot * nu(1)   0];

% System mass and Coriolis-centripetal matrices
M = MRB + MA;
C = CRB + CA;

% Linear damping terms (hydrodynamic derivatives)
Xu = -24.4 * g / Umax;        % specified using the maximum speed  
Yv = -M(2,2) / T_sway;        % specified using the time constant in sway
Nr = -M(3,3) / T_yaw;         % specified using the time constant in T_yaw

% 2-DOF constant input matrix B_prop for the propellers in sway and yaw
B_prop = k_pos * [1, 1; y_pont, -y_pont];

% Control forces and moments, with saturated propeller speed
Thrust = zeros(2,1);
for i = 1:1:2                  
     Thrust(i) = n(i);    %  thrust (N) 
end


% Control forces and moments
tau = [Thrust(1) + Thrust(2); 0; -l1 * Thrust(1) - l2 * Thrust(2)];

% Linear damping using relative velocities + nonlinear yaw damping
Xh = Xu * nu(1);
Yh = Yv * nu(2); 
Nh = Nr * (1 + 10 * abs(nu(3))) * nu(3);

tau_damp = [Xh; Yh; Nh];

% Trim condition: G * eta_0 = g_0
eta_0 = [0; 0; 0];
eta = eta - eta_0; % shifted equilibrium

% Kinematic transformation matrix
J = [cos(eta(3)), -sin(eta(3)), 0;
     sin(eta(3)),  cos(eta(3)), 0;
     0,            0,           1];

% Time derivative of the state vector
xdot = [J * nu;                 % derivative of position [xdot, ydot, psidot]
        M \ (tau + tau_damp - C * nu)];  % derivative of velocities [udot, vdot, rdot]
     
end