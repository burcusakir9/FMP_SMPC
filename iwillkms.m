x_eq = [5; 2; 0; 0; 0; 0];
u_eq = [0; 0];
Ts = 0.1;

waypoint = [10; 10];

funnel.rMax = 8.0;
funnel.rMin = 0.0;
funnel.betaMax = deg2rad(20);

limits.uMin  = [-8; -8];
limits.uMax  = [ 8;  8];
limits.duMin = [-1; -1];
limits.duMax = [ 1;  1];

limits.xMin = [-inf; -inf; -inf; -2; -2; -1];
limits.xMax = [ inf;  inf;  inf;  2;  2;  1];

mpcOptions.constraintMode = 'polar_linearized';

mdl = tugboat3d_linear_model(x_eq, u_eq, Ts, waypoint, funnel, limits, mpcOptions);
disp(mdl.A)
disp(mdl.B)
disp(mdl.Hx)
disp(mdl.hx)

function [Hf, hf, info] = build_polar_funnel_constraint(x_eq, waypoint, funnel)
% Polar-coordinate linearized constraints around x_eq
%
% rho = sqrt((X-Xw)^2 + (Y-Yw)^2)
% theta = atan2(Y-Yw, X-Xw)
%
% Constraints example:
%   rho <= rMax
%   rho >= rMin        (optional)
%   |theta - theta_eq| <= betaMax   (optional local angular funnel)
%
% All are linearized in deviation coordinates dx.

    Xeq = x_eq(1);
    Yeq = x_eq(2);

    Xw = waypoint(1);
    Yw = waypoint(2);

    ex = Xeq - Xw;
    ey = Yeq - Yw;

    rho_eq = sqrt(ex^2 + ey^2);

    if rho_eq < 1e-8
        error('Polar linearization is singular because x_eq is at the waypoint.');
    end

    if ~isfield(funnel, 'rMin') || isempty(funnel.rMin)
        funnel.rMin = 0;
    end

    Hf = [];
    hf = [];

    %--------------------------------------------------
    % 1) rho <= rMax
    % rho linearization:
    %   rho ≈ rho_eq + drho_dX*dX + drho_dY*dY
    % drho = [ex/rho_eq, ey/rho_eq]
    % therefore:
    %   drho_grad * [dxX; dxY] <= rMax - rho_eq
    %--------------------------------------------------
    drho = [ex/rho_eq; ey/rho_eq];

    Hrho_max = zeros(1,6);
    Hrho_max(1,1) = drho(1);
    Hrho_max(1,2) = drho(2);

    hrho_max = funnel.rMax - rho_eq;

    Hf = [Hf; Hrho_max];
    hf = [hf; hrho_max];

    %--------------------------------------------------
    % 2) rho >= rMin   ->  -rho <= -rMin
    %--------------------------------------------------
    if funnel.rMin > 0
        Hrho_min = zeros(1,6);
        Hrho_min(1,1) = -drho(1);
        Hrho_min(1,2) = -drho(2);

        hrho_min = -(funnel.rMin - rho_eq);

        Hf = [Hf; Hrho_min];
        hf = [hf; hrho_min];
    end

    %--------------------------------------------------
    % 3) Angular funnel:
    %    |theta - theta_eq| <= betaMax
    %
    % theta = atan2(ey, ex)
    % dtheta/dX = -ey/(ex^2+ey^2)
    % dtheta/dY =  ex/(ex^2+ey^2)
    %
    % local linearized form:
    %    theta ≈ theta_eq + grad_theta * [dxX; dxY]
    %
    % so:
    %    grad_theta * [dxX; dxY] <= betaMax
    %   -grad_theta * [dxX; dxY] <= betaMax
    %--------------------------------------------------
    theta_eq = atan2(ey, ex);

    if isfield(funnel, 'betaMax') && ~isempty(funnel.betaMax)
        gtheta = [-ey/(rho_eq^2); ex/(rho_eq^2)];

        Hth1 = zeros(1,6);
        Hth1(1,1) =  gtheta(1);
        Hth1(1,2) =  gtheta(2);

        Hth2 = zeros(1,6);
        Hth2(1,1) = -gtheta(1);
        Hth2(1,2) = -gtheta(2);

        hth1 = funnel.betaMax;
        hth2 = funnel.betaMax;

        Hf = [Hf; Hth1; Hth2];
        hf = [hf; hth1; hth2];
    else
        gtheta = [];
    end

    info.type = 'polar_linearized';
    info.rho_eq = rho_eq;
    info.theta_eq = theta_eq;
    info.drho = drho;
    info.gtheta = gtheta;
    info.note = 'Local linearization of polar geometry around equilibrium';
end

function [Hf, hf, info] = build_cartesian_funnel_constraint(x_eq, waypoint, funnel)
% Linearize circular constraint around x_eq
%
% Exact outer circle:
%   c(X,Y) = (X-Xw)^2 + (Y-Yw)^2 - rMax^2 <= 0
%
% Linearized at equilibrium:
%   c_eq + grad_c' * [dX; dY] <= 0
%
% Since dx = x - x_eq:
%   grad_c' * [dx_X; dx_Y] <= -c_eq
%
% This gives a local linear inequality.

    Xeq = x_eq(1);
    Yeq = x_eq(2);

    Xw = waypoint(1);
    Yw = waypoint(2);

    rMax = funnel.rMax;

    ex = Xeq - Xw;
    ey = Yeq - Yw;

    c_eq = ex^2 + ey^2 - rMax^2;
    grad = [2*ex; 2*ey];

    Hf = zeros(1,6);
    Hf(1,1) = grad(1);
    Hf(1,2) = grad(2);

    hf = -c_eq;

    info.type = 'cartesian_linearized_outer_circle';
    info.c_eq = c_eq;
    info.grad = grad;
    info.note = 'Local tangent approximation of outer circle';
end

function [Hx, hx] = build_box_state_constraints(x_eq, limits)
% Build deviation-state constraints:
%   xMin <= x_eq + dx <= xMax

    Hx = [];
    hx = [];

    hasMin = isfield(limits, 'xMin') && ~isempty(limits.xMin);
    hasMax = isfield(limits, 'xMax') && ~isempty(limits.xMax);

    if ~hasMin && ~hasMax
        return;
    end

    n = length(x_eq);
    Hx = [];
    hx = [];

    if hasMax
        Hx = [Hx; eye(n)];
        hx = [hx; limits.xMax - x_eq];
    end

    if hasMin
        Hx = [Hx; -eye(n)];
        hx = [hx; -(limits.xMin - x_eq)];
    end
end

function [Hdu, hdu] = build_input_rate_constraints(limits)
% Build:
%   duMin <= Delta u <= duMax
% => [ I; -I ] Delta u <= [duMax; -duMin]

    if ~isfield(limits, 'duMin') || ~isfield(limits, 'duMax')
        Hdu = [];
        hdu = [];
        return;
    end

    I = eye(2);

    Hdu = [ I;
           -I];

    hdu = [limits.duMax;
          -limits.duMin];
end

function [Hu, hu] = build_input_constraints(u_eq, limits)
% Build constraints on deviation input du:
%   uMin <= u_eq + du <= uMax
% => [ I; -I ] du <= [uMax-u_eq; -(uMin-u_eq)]

    if ~isfield(limits, 'uMin') || ~isfield(limits, 'uMax')
        Hu = [];
        hu = [];
        return;
    end

    I = eye(2);

    Hu = [ I;
          -I];

    hu = [limits.uMax - u_eq;
         -(limits.uMin - u_eq)];
end

function [A, B] = tugboat3d_linearized_jacobian(x_eq, u_eq)
% Strict Jacobian linearization around (x_eq, u_eq)
%
% State:
%   x = [X; Y; psi; ub; vb; r]
% Input:
%   u = [FL; FR]

    psi = x_eq(3);
    ub  = x_eq(4);
    vb  = x_eq(5);
    r   = x_eq(6);

    %#ok<NASGU>
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

    M = [m - Xudot,      0,            0;
         0,              m - Yvdot,   -Yrdot;
         0,             -Yrdot,        Iz - Nrdot];

    m11 = M(1,1);
    m22 = M(2,2);
    m23 = M(2,3);

    A = zeros(6,6);

    cpsi = cos(psi);
    spsi = sin(psi);

    % Kinematics
    A(1,3) = -spsi*ub - cpsi*vb;
    A(1,4) =  cpsi;
    A(1,5) = -spsi;

    A(2,3) =  cpsi*ub - spsi*vb;
    A(2,4) =  spsi;
    A(2,5) =  cpsi;

    A(3,6) = 1;

    % Dynamics Jacobian dF/dnu
    dF_dnu = zeros(3,3);

    % d(|x|x)/dx = 2|x| for x ~= 0
    % at x = 0 this is not differentiable; use 0 as local approximation
    dub = 2*abs(ub);
    dvb = 2*abs(vb);
    dr  = 2*abs(r);

    % Row 1
    dF_dnu(1,1) = Xu + Xuu*dub;
    dF_dnu(1,2) = m22*r;
    dF_dnu(1,3) = m22*vb + 2*m23*r;

    % Row 2
    dF_dnu(2,1) = -m11*r;
    dF_dnu(2,2) = Yv + Yvv*dvb;
    dF_dnu(2,3) = Yr - m11*ub;

    % Row 3
    dF_dnu(3,1) = -(m22 - m11)*vb - m23*r;
    dF_dnu(3,2) = Nv - (m22 - m11)*ub;
    dF_dnu(3,3) = Nr + Nrr*dr - m23*ub;

    A_nu = M \ dF_dnu;
    A(4:6,4:6) = A_nu;

    % Input Jacobian
    B_prop = [1, 1;
              0, 0;
              d, -d];

    B = zeros(6,2);
    B(4:6,:) = M \ B_prop;
end

function xdot = tugboat3d_dynamics(x, u)
% Nonlinear 3-DOF tugboat model
%
% x = [X; Y; psi; ub; vb; r]
% u = [FL; FR]

    psi = x(3);
    ub  = x(4);
    vb  = x(5);
    r   = x(6);

    FL = u(1);
    FR = u(2);

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

    % Kinematics
    Xdot   = cos(psi)*ub - sin(psi)*vb;
    Ydot   = sin(psi)*ub + cos(psi)*vb;
    psidot = r;

    % Dynamics
    F1 = FL + FR ...
       + Xu*ub + Xuu*abs(ub)*ub ...
       + m22*vb*r + m23*r^2;

    F2 = Yv*vb + Yr*r + Yvv*abs(vb)*vb ...
       - m11*ub*r;

    F3 = d*(FL - FR) ...
       + Nv*vb + Nr*r + Nrr*abs(r)*r ...
       - (m22 - m11)*ub*vb - m23*ub*r;

    nudot = M \ [F1; F2; F3];

    xdot = [Xdot; Ydot; psidot; nudot];
end


function mdl = tugboat3d_linear_model(x_eq, u_eq, Ts, waypoint, funnel, limits, options)
% Strictly linear MPC model for tugboat
%
% Linearization is around an equilibrium/trim point:
%       f(x_eq, u_eq) = 0
%
% Model is in deviation variables:
%       dx = x - x_eq
%       du = u - u_eq
%
% Continuous-time:
%       d/dt(dx) = A*dx + B*du
%
% Discrete-time:
%       dx(k+1) = Ad*dx(k) + Bd*du(k)
%
% INPUTS:
%   x_eq     : [6x1] equilibrium state
%   u_eq     : [2x1] equilibrium input
%   Ts       : sample time
%   waypoint : [Xw; Yw]
%   funnel   : struct with fields:
%                .rMax
%                .rMin      (optional, default 0)
%                .betaMax   (optional, only for polar option)
%   limits   : struct with fields:
%                .uMin, .uMax
%                .duMin, .duMax
%                .xMin, .xMax   (optional)
%   options  : struct with field:
%                .constraintMode = 'cartesian_linearized' or 'polar_linearized'
%
% OUTPUT:
%   mdl.A, mdl.B          continuous-time strict linear model
%   mdl.Ad, mdl.Bd        discrete-time strict linear model
%   mdl.Hx, mdl.hx        state constraints in deviation variables: Hx*dx <= hx
%   mdl.Hu, mdl.hu        input constraints in deviation variables: Hu*du <= hu
%   mdl.Hdu, mdl.hdu      input-rate constraints: Hdu*ddu <= hdu
%   mdl.Aaug, mdl.Baug    augmented model for rate-constrained MPC
%   mdl.x_eq, mdl.u_eq
%   mdl.info

    arguments
        x_eq (6,1) double
        u_eq (2,1) double
        Ts (1,1) double {mustBePositive}
        waypoint (2,1) double
        funnel struct
        limits struct
        options.constraintMode char {mustBeMember(options.constraintMode,...
            {'cartesian_linearized','polar_linearized'})} = 'cartesian_linearized'
    end

    % Check equilibrium consistency
    f_eq = tugboat3d_dynamics(x_eq, u_eq);
    if norm(f_eq, 2) > 1e-6
        warning('x_eq, u_eq are not an exact equilibrium. Norm(f_eq) = %.3e', norm(f_eq,2));
        warning('The resulting model is mathematically a local linearization in deviation variables,');
        warning('but for a strictly linear offset-free model you should use a true trim point.');
    end

    % Linearization
    [A, B] = tugboat3d_linearized_jacobian(x_eq, u_eq);

    % Discretization (Euler)
    % For better accuracy you may later replace this by c2d
    Ad = eye(6) + Ts*A;
    Bd = Ts*B;

    % Input constraints in deviation variables:
    % uMin <= u_eq + du <= uMax
    [Hu, hu] = build_input_constraints(u_eq, limits);

    % Input-rate constraints:
    % duMin <= Delta u <= duMax
    [Hdu, hdu] = build_input_rate_constraints(limits);

    % State constraints in deviation variables:
    % xMin <= x_eq + dx <= xMax
    [Hx_box, hx_box] = build_box_state_constraints(x_eq, limits);

    % Funnel constraints
    switch options.constraintMode
        case 'cartesian_linearized'
            [Hf, hf, infoF] = build_cartesian_funnel_constraint(x_eq, waypoint, funnel);

        case 'polar_linearized'
            [Hf, hf, infoF] = build_polar_funnel_constraint(x_eq, waypoint, funnel);
    end

    % Combine state constraints
    Hx = [Hx_box; Hf];
    hx = [hx_box; hf];

    % Augmented model for input-rate-constrained MPC
    %
    % Let z = [dx; du_prev]
    % Let v = Delta u = du(k) - du(k-1)
    %
    % Then:
    %   dx(k+1)      = Ad dx(k) + Bd du(k)
    %                = Ad dx(k) + Bd(du_prev(k) + v(k))
    %
    %   du_prev(k+1) = du_prev(k) + v(k)
    %
    % Therefore:
    %   z(k+1) = Aaug z(k) + Baug v(k)
    %
    Aaug = [Ad, Bd;
            zeros(2,6), eye(2)];

    Baug = [Bd;
            eye(2)];

    mdl = struct();
    mdl.A = A;
    mdl.B = B;
    mdl.Ad = Ad;
    mdl.Bd = Bd;

    mdl.Hx = Hx;
    mdl.hx = hx;

    mdl.Hu = Hu;
    mdl.hu = hu;

    mdl.Hdu = Hdu;
    mdl.hdu = hdu;

    mdl.Aaug = Aaug;
    mdl.Baug = Baug;

    mdl.x_eq = x_eq;
    mdl.u_eq = u_eq;

    mdl.info.f_eq = f_eq;
    mdl.info.funnel = infoF;
    mdl.info.constraintMode = options.constraintMode;
end