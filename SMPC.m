%% SMPC (Stochastic MPC) for the thesis-style 2D holonomic robot
% - Linear discrete model: [px py vx vy]'
% - Inputs: [ax ay]' (accelerations)
% - Process noise: w ~ N(0,Q)
% - Measurement: y = [px py] + v, v ~ N(0,R)
% - KF for state estimation
% - SMPC via chance-constraint tightening using predicted covariance
%
% This is a self-contained script. Replace the "waypoints" with your node/overlap points later.

clear; clc; close all;
rng(7);

%% ------------------- Parameters -------------------
dt = 0.05;                 % [s]
N  = 10;                  % MPC horizon
T  = 400;                 % simulation steps (one run)
MC = 1;                  % Monte Carlo runs

% Robot model (double integrator in x,y)
A = [1 0 dt 0;
     0 1 0  dt;
     0 0 1  0;
     0 0 0  1];

B = [0.5*dt^2 0;
     0 0.5*dt^2;
     dt 0;
     0 dt];

C = [1 0 0 0;
     0 1 0 0];

nx = size(A,1);
nu = size(B,2);
ny = size(C,1);

% Noise (tune to match thesis-like levels)
Q = diag([1e-6, 1e-6, 1e-6, 1e-6]);   % process noise covariance
R = diag([1e-3, 1e-3]);       % measurement noise covariance

% Precompute Cholesky factors for noise generation
L_Q = chol(Q, 'lower'); 
L_R = chol(R, 'lower');

% KF initial covariance
P0 = diag([1e-2 1e-2 5e-2 5e-2]);

% Constraints (example)
xmin = 0; xmax = 16;
ymin = 0; ymax = 8;

umax = 3.0;      % max accel per axis (hard limit in deterministic MPC)

% Chance constraints (tighten mean constraints)
beta_state = 0.05;   % probability of violating state bounds (per step, per axis, two-sided)
beta_input = 0.005;   % probability of violating input bounds (per step, per axis, two-sided)

% Chance constraints (tighten mean constraints) using erfinv (Core MATLAB)
% k_state = norminv(1 - beta_state/2);  
p_state = 1 - beta_state/2; % constraint is split into two for x and y
k_state = sqrt(2) * erfinv(2 * p_state - 1);  %vinverse CDF

% k_input = norminv(1 - beta_input/2);  
p_input = 1 - beta_input/2;
k_input = sqrt(2) * erfinv(2 * p_input - 1);

% Cost weights
Qpos = 100;    % position tracking
Qvel = 2;     % velocity penalty
Ru   = 0.1;   % input penalty
W = diag([Qpos Qpos Qvel Qvel]);   % stage weight on x
Wu = Ru*eye(nu);                   % stage weight on u
WN = 10*W;                         % terminal weight

% Waypoints (replace with your node-intersection centroid points later)
waypoints = [ 1.0  1.0;
              2.0  6.0;
              6.0  6.0;
              9.5  2.0;
              13.5 6.5;
              15.0 5.0 ];
wp_idx = 1;
wp_tol = 0.25;   % switch waypoint when within this distance

% Initial true state and initial estimate
x_true0 = [waypoints(1,1); waypoints(1,2); 0; 0];
xhat0   = x_true0 + [0.05; -0.05; 0.1; -0.1];

%% ------------------- Precompute MPC prediction matrices -------------------
% X = Sx*x0 + Su*U, where X stacks x1..xN (or x0..xN depending convention).
[Sx, Su] = predMatrices(A,B,N);

% Build block-diagonal weights for stacked horizon:
% We optimize U for cost sum_{i=0..N-1} x_i'W x_i + u_i'Wu u_i + x_N'WN x_N
% Here we stack states x1..xN in X, and use x0 separately.
Wx = blkdiag(kron(eye(N-1), W), WN);     % for x1..xN (N blocks, last is WN)
Wu_blk = kron(eye(N), Wu);              % for u0..u_{N-1}

% QP form: 0.5*U'HU + f'U
H = Su'*Wx*Su + Wu_blk;      % symmetric PSD
H = (H+H')/2;                % numerical sym

optsQP = optimoptions('quadprog','Display','off');

%% ------------------- Monte Carlo simulation -------------------
traj_all = cell(MC,1);

for mc = 1:MC
    x_true = x_true0;
    xhat = xhat0;
    P = P0;

    wp_idx = 1;
    U_hist = zeros(nu,T);
    X_hist = zeros(nx,T);
    Xh_hist = zeros(nx,T);

    for k = 1:T
        % Choose current waypoint reference
        r = waypoints(wp_idx,:)';
        if norm(x_true(1:2) - r) < wp_tol && wp_idx < size(waypoints,1)
            wp_idx = wp_idx + 1;
            r = waypoints(wp_idx,:)';
        end

        % Build stacked reference trajectory (constant reference over horizon)
        Xref = repmat([r; 0; 0], N, 1);   % desire stop at waypoint (vx=vy=0)

        % Predict covariance over horizon for tightening
        % P_i for i=1..N (matching X = x1..xN)
        Ppred = zeros(nx,nx,N);
        Ptmp = P;
        for i = 1:N
            Ptmp = A*Ptmp*A' + Q;
            Ppred(:,:,i) = Ptmp;
        end

        % Build tightened linear constraints for state bounds on mean:
        % For each i, enforce:
        %   xmin <= mu_px_i <= xmax, ymin <= mu_py_i <= ymax with tightening
        % mu_i = (Sx*xhat + Su*U)_i
        %
        % Tightening per bound using std dev from Ppred:
        % mu_px_i >= xmin + k_state*sqrt(var_px_i)
        % mu_px_i <= xmax - k_state*sqrt(var_px_i)
        %
        % Same for y. (No constraints on v here; easy to add.)
        Aineq = [];
        bineq = [];

        % Extract rows that map U -> px_i, py_i in stacked X
        % X stack order: [x1; x2; ...; xN], each x is 4x1
        for i = 1:N
            row_px = (i-1)*nx + 1;
            row_py = (i-1)*nx + 2;

            Su_px = Su(row_px, :);
            Su_py = Su(row_py, :);

            mu_px_aff = Sx(row_px,:)*xhat;
            mu_py_aff = Sx(row_py,:)*xhat;

            sig_px = sqrt(max(Ppred(1,1,i),0));
            sig_py = sqrt(max(Ppred(2,2,i),0));

            lb_px = xmin + k_state*sig_px;
            ub_px = xmax - k_state*sig_px;
            lb_py = ymin + k_state*sig_py;
            ub_py = ymax - k_state*sig_py;

            % mu_px = mu_px_aff + Su_px*U
            % enforce mu_px <= ub_px  => Su_px*U <= ub_px - mu_px_aff
            Aineq = [Aineq;  Su_px];
            bineq = [bineq;  ub_px - mu_px_aff];

            % enforce -mu_px <= -lb_px => -Su_px*U <= -(lb_px - mu_px_aff)
            Aineq = [Aineq; -Su_px];
            bineq = [bineq; -(lb_px - mu_px_aff)];

            % mu_py
            Aineq = [Aineq;  Su_py];
            bineq = [bineq;  ub_py - mu_py_aff];

            Aineq = [Aineq; -Su_py];
            bineq = [bineq; -(lb_py - mu_py_aff)];
        end

        % Tightened input bounds on mean:
        % |u_j| <= umax - k_input*sig_u
        % Here we approximate sig_u = 0 (no input noise model) -> deterministic.
        % If you have input noise, add it similarly.
        U_lb = -umax*ones(nu*N,1);
        U_ub =  umax*ones(nu*N,1);

        % Linear term f = Su'Wx(Sx*xhat - Xref)
        f = 2 * Su' * Wx * (Sx * xhat - Xref);

        % Solve QP
        [Uopt,~,exitflag] = quadprog(H, f, Aineq, bineq, [], [], U_lb, U_ub, [], optsQP);

        if exitflag <= 0 || isempty(Uopt)
            % If infeasible, apply zero input (or last input)
            u = zeros(nu,1);
        else
            u = Uopt(1:nu);
        end

        % Apply to true system with process noise
        % Apply to true system with process noise (Replaces mvnrnd)
        w = L_Q * randn(nx, 1);
        x_true = A*x_true + B*u + w;
        
        % Measurement (Replaces mvnrnd)
        v = L_R * randn(ny, 1);
        y = C*x_true + v;

        % Kalman Filter update
        % Predict
        xpred = A*xhat + B*u;
        PpredKF = A*P*A' + Q;
        % Update
        S = C*PpredKF*C' + R;
        K = PpredKF*C' / S;
        xhat = xpred + K*(y - C*xpred);
        P = (eye(nx) - K*C)*PpredKF;

        % log
        U_hist(:,k) = u;
        X_hist(:,k) = x_true;
        Xh_hist(:,k)= xhat;

        % stop if final waypoint reached
        if wp_idx == size(waypoints,1) && norm(x_true(1:2) - waypoints(end,:)') < wp_tol
            % fill rest for plotting consistency
            X_hist(:,k+1:end) = repmat(x_true,1,T-k);
            Xh_hist(:,k+1:end)= repmat(xhat,1,T-k);
            break;
        end
    end

    traj_all{mc} = X_hist;
end

%% ------------------- Plot -------------------
figure('Color','w'); hold on; axis equal; grid on;
xlabel('x'); ylabel('y'); title('SMPC Monte Carlo trajectories (holonomic robot)');

% plot workspace bounds
plot([xmin xmax xmax xmin xmin],[ymin ymin ymax ymax ymin],'k-','LineWidth',1);

% plot waypoints
plot(waypoints(:,1), waypoints(:,2), 'ko','LineWidth',1.5);

% plot trajectories
for mc = 1:MC
    X = traj_all{mc};
    plot(X(1,:), X(2,:));
end

legend({'Workspace','Waypoints'}, 'Location','best');

%% ------------------- Helpers -------------------
function [Sx, Su] = predMatrices(A,B,N)
% Build prediction matrices for X = [x1; x2; ...; xN] = Sx*x0 + Su*U
% U = [u0; u1; ...; u_{N-1}]
    nx = size(A,1);
    nu = size(B,2);

    Sx = zeros(nx*N, nx);
    Su = zeros(nx*N, nu*N);

    A_pow = eye(nx);
    for i = 1:N
        A_pow = A_pow * A;                 % A^i
        Sx((i-1)*nx+1:i*nx, :) = A_pow;
        for j = 1:i
            % contribution of u_{j-1} to x_i is A^{i-j} B
            A_ij = A^(i-j);
            Su((i-1)*nx+1:i*nx, (j-1)*nu+1:j*nu) = A_ij * B;
        end
    end
end
