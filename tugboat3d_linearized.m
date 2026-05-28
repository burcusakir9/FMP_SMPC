function [A, B] = tugboat3d_linearized(x_eq, u_eq)
% Analytical Jacobian of tugboat3d around (x_eq, u_eq)
%
% x = [X; Y; psi; ub; vb; r]
% u = [FL; FR]

    X   = x_eq(1);
    Y   = x_eq(2);
    psi = x_eq(3);
    ub  = x_eq(4);
    vb  = x_eq(5);
    r   = x_eq(6);

    FL = u_eq(1);
    FR = u_eq(2);

    % Parameters
    m  = 10.2;
    Iz = 0.63994;

    beam = 0.29;
    d = beam/2;

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

    % Inertia matrix
    M = [m - Xudot,      0,            0;
         0,              m - Yvdot,   -Yrdot;
         0,             -Yrdot,        Iz - Nrdot];

    m11 = M(1,1);
    m22 = M(2,2);
    m23 = M(2,3);

    % ---------- Kinematics Jacobian ----------
    cpsi = cos(psi);
    spsi = sin(psi);

    A = zeros(6,6);

    % Xdot = cos(psi)*ub - sin(psi)*vb
    A(1,3) = -spsi*ub - cpsi*vb;
    A(1,4) =  cpsi;
    A(1,5) = -spsi;

    % Ydot = sin(psi)*ub + cos(psi)*vb
    A(2,3) =  cpsi*ub - spsi*vb;
    A(2,4) =  spsi;
    A(2,5) =  cpsi;

    % psidot = r
    A(3,6) = 1;

    % ---------- Dynamics Jacobian ----------
    % F = tau + tau_h - C*nu
    % nudot = M^{-1} F
    %
    % With:
    % F1 = FL+FR + Xu*ub + Xuu*|ub|ub + m22*vb*r + m23*r^2
    % F2 =         Yv*vb + Yr*r + Yvv*|vb|vb - m11*ub*r
    % F3 = d(FL-FR)+ Nv*vb + Nr*r + Nrr*|r|r - (m22-m11)*ub*vb - m23*ub*r

    dabs_ub = 2*abs(ub);
    dabs_vb = 2*abs(vb);
    dabs_r  = 2*abs(r);

    dF_dnu = zeros(3,3);

    % Row 1 derivatives
    dF_dnu(1,1) = Xu + Xuu*dabs_ub;
    dF_dnu(1,2) = m22*r;
    dF_dnu(1,3) = m22*vb + 2*m23*r;

    % Row 2 derivatives
    dF_dnu(2,1) = -m11*r;
    dF_dnu(2,2) = Yv + Yvv*dabs_vb;
    dF_dnu(2,3) = Yr - m11*ub;

    % Row 3 derivatives
    dF_dnu(3,1) = -(m22 - m11)*vb - m23*r;
    dF_dnu(3,2) = Nv - (m22 - m11)*ub;
    dF_dnu(3,3) = Nr + Nrr*dabs_r - m23*ub;

    A_nu = M \ dF_dnu;

    A(4:6,4:6) = A_nu;

    % ---------- Input Jacobian ----------
    B_prop = [1, 1;
              0, 0;
              d, -d];

    B = zeros(6,2);
    B(4:6,:) = M \ B_prop;
end