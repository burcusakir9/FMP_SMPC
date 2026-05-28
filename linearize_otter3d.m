function [A, B] = linearize_otter3d(x_eq, u_eq)
    % Linearize the Otter USV 6-state model around (x_eq, u_eq)
    % Consistent with state order [x, y, yaw, u, v, r]
    
    % Extract parameters from otter3d_reordered
    [~, ~, M, B_prop] = otter3d(zeros(6,1), zeros(2,1));
    
    % Main parameters
    m = M(1,1);         % Mass
    Iz = M(3,3);        % Moment of inertia
    Xudot = -M(1,1);    % Added mass in surge
    Yvdot = -M(2,2);    % Added mass in sway
    Nrdot = -M(3,3);    % Added mass in yaw
    
    % From otter3d_reordered
    g = 9.81;
    Umax = 6 * 0.5144;
    Xu = -24.4 * g / Umax;
    Yv = -M(2,2);       % T_sway = 1
    Nr = -M(3,3);       % T_yaw = 1
    k = 0.02216/2;      % Thruster coefficient
    l = 0.395;          % Lever arm (y_pont)
    % Symbolic variables - note order is [u, v, r] for velocities
    syms u v r real
    syms u1 u2 real
    nu = [u; v; r]; 
    % --- Inertia matrix ---
    MRB = diag([m, m, Iz]);
    MA = -diag([Xudot, Yvdot, Nrdot]);
    M = MRB + MA;
    % --- Coriolis matrices ---
    CRB = [  0,   0,  -m*v;
             0,   0,   m*u;
             m*v, -m*u, 0 ];
    CA = [  0,      0,      Yvdot*v;
            0,      0,     -Xudot*u;
           -Yvdot*v, Xudot*u,  0 ];
    C = CRB + CA;
    % --- Damping matrix ---
    D = [ Xu,   0,    0;
           0,   Yv,   0;
           0,    0,   Nr*(1 + 10*abs(r)) ];
    % --- Dynamics f(nu) = M⁻¹ (D*nu - C*nu) ---
    f_nu = M \ (D*nu - C*nu);
    
    % Linearize about equilibrium velocities (x_eq(4:6) = [u, v, r]
    vel_eq = x_eq(4:6);
    A_sym = jacobian(f_nu, [u; v; r]);
    A_vel = double(subs(A_sym, [u; v; r], vel_eq));
    % For position states: [x, y, yaw] dynamics are linear
    J = [cos(x_eq(3)), -sin(x_eq(3)), 0;
         sin(x_eq(3)),  cos(x_eq(3)), 0;
         0,            0,            1];
    
    % Construct full A matrix
    A = [zeros(3,3),    J;
         zeros(3,3), A_vel];
    % --- Input matrix ---
    % Thruster forces (same as otter3d_reordered)
    tau = [ k*(u1*abs(u1) + u2*abs(u2));
            0;
            l*k*(u1*abs(u1) - u2*abs(u2)) ];
    
    g_u = M \ tau;
    B_sym = jacobian(g_u, [u1; u2]);
    B_vel = double(subs(B_sym, [u1; u2], u_eq));
    % Position states are not directly actuated
    B = [zeros(3,2);
         B_vel];
end 