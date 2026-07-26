function [A_c, B_c] = tugboat3d_linearized(x_bar, u_bar)
%TUGBOAT3D_LINEARIZED Continuous-time analytical linearization.
%
% Linearizes tugboat3d around the nominal state-input pair:
%
%   x_bar = [X; Y; psi; ub; vb; r]
%   u_bar = [FL; FR]
%
% Returns the continuous-time error-state model:
%
%   delta_xdot = A_c*delta_x + B_c*delta_u
%
% where:
%
%   delta_x = x - x_bar
%   delta_u = u - u_bar
%
% The matrices can subsequently be discretized for MPC.

    % Ensure column vectors
    x_bar = x_bar(:);
    u_bar = u_bar(:);

    % Validate dimensions
    if numel(x_bar) ~= 6
        error('x_bar must contain 6 elements.');
    end

    if numel(u_bar) ~= 2
        error('u_bar must contain 2 elements.');
    end

    % State at the linearization point
    psi = x_bar(3);
    ub  = x_bar(4);
    vb  = x_bar(5);
    r   = x_bar(6);

    % u_bar does not currently affect A_c or B_c because the nonlinear
    % model is affine in FL and FR. It is retained in the interface so that
    % the linearization point is represented consistently.
    %#ok<NASGU>

    %%%%%%%%%%%%%%%%%%%%%%%
    % Model parameters    %
    %%%%%%%%%%%%%%%%%%%%%%%

    m  = 10.2;
    Iz = 0.63994;

    beam = 0.29;

    % This assumes that each thrust line is located at beam/2 from the
    % vessel centerline. Replace with the actual CG-to-thrust-line distance
    % when available.
    d = beam / 2.0;

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

    %%%%%%%%%%%%%%%%%%%%%%%
    % Inertia matrix      %
    %%%%%%%%%%%%%%%%%%%%%%%

    % Assumes added-mass symmetry:
    % Nvdot = Yrdot
    M = [m - Xudot,  0.0,          0.0;
         0.0,        m - Yvdot,   -Yrdot;
         0.0,       -Yrdot,        Iz - Nrdot];

    m11 = M(1,1);
    m22 = M(2,2);
    m23 = M(2,3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Continuous-time state Jacobian    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A_c = zeros(6,6);

    cpsi = cos(psi);
    spsi = sin(psi);

    % Xdot = cos(psi)*ub - sin(psi)*vb
    A_c(1,3) = -spsi*ub - cpsi*vb;
    A_c(1,4) =  cpsi;
    A_c(1,5) = -spsi;

    % Ydot = sin(psi)*ub + cos(psi)*vb
    A_c(2,3) =  cpsi*ub - spsi*vb;
    A_c(2,4) =  spsi;
    A_c(2,5) =  cpsi;

    % psidot = r
    A_c(3,6) = 1.0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dynamic part of state Jacobian    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % nudot = M^(-1) F(nu,u)
    %
    % F1 = FL + FR
    %      + Xu*ub + Xuu*abs(ub)*ub
    %      + m22*vb*r + m23*r^2
    %
    % F2 = Yv*vb + Yr*r + Yvv*abs(vb)*vb
    %      - m11*ub*r
    %
    % F3 = d*(FL-FR)
    %      + Nv*vb + Nr*r + Nrr*abs(r)*r
    %      - (m22-m11)*ub*vb
    %      - m23*ub*r

    % For g(q) = abs(q)*q:
    %
    % dg/dq = 2*abs(q)
    dub_drag = 2.0 * abs(ub);
    dvb_drag = 2.0 * abs(vb);
    dr_drag  = 2.0 * abs(r);

    dF_dnu = zeros(3,3);

    % F1 derivatives with respect to [ub, vb, r]
    dF_dnu(1,1) = Xu + Xuu*dub_drag;
    dF_dnu(1,2) = m22*r;
    dF_dnu(1,3) = m22*vb + 2.0*m23*r;

    % F2 derivatives with respect to [ub, vb, r]
    dF_dnu(2,1) = -m11*r;
    dF_dnu(2,2) = Yv + Yvv*dvb_drag;
    dF_dnu(2,3) = Yr - m11*ub;

    % F3 derivatives with respect to [ub, vb, r]
    dF_dnu(3,1) = -(m22 - m11)*vb - m23*r;
    dF_dnu(3,2) = Nv - (m22 - m11)*ub;
    dF_dnu(3,3) = Nr + Nrr*dr_drag - m23*ub;

    A_c(4:6,4:6) = M \ dF_dnu;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Continuous-time input Jacobian    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    B_prop = [1.0,  1.0;
              0.0,  0.0;
              d,   -d];

    B_c = zeros(6,2);
    B_c(4:6,:) = M \ B_prop;
end