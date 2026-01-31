%% SMPC (Stochastic MPC) for the thesis-style 2D holonomic robot
% - Linear discrete model: [px py vx vy]'
% - Inputs: [ax ay]' (accelerations)
% - Process noise: w ~ N(0,Q)
% - MeaHrement: y = [px py] + v, v ~ N(0,R)
% - KF for state estimation
% - SMPC via chance-constraint tightening using predicted covariance

clear; clc; close all;
rng(7);

%% ------------------- Parameters -------------------
dt = 0.05; %sampling time, seconds
N  = 10;  % MPC horizon
T  = 400; % simulation steps 
MC = 1; % Monte Carlo runs

% Robot model
A = [1 0 dt 0;
     0 1 0  dt;
     0 0 1  0;
     0 0 0  1];

B = [0.5*dt^2 0;
     0 0.5*dt^2;
     dt 0;
     0 dt];

C = [1 0 0 0;
     0 1 0 0]; % position states are observable

nx = size(A,1);
nu = size(B,2);
ny = size(C,1);

Q = diag([1e-6, 1e-6, 1e-6, 1e-6]);   % process noise covariance
R = diag([1e-6, 1e-6]);       % meaHrement noise covariance
Qu = diag([0.05, 0.05]); % Input noise

% Compute Cholesky factors in order to add noise,  Q = L_Q * L_Q'
L_Q = chol(Q, 'lower'); 
L_R = chol(R, 'lower');
L_Qu = chol(Qu, 'lower');

% KF initial covariance
P0 = diag([1e-6 1e-6 1e-6 1e-6]);

% Map limits, linear constraints
xmin = 0; xmax = 16;
ymin = 0; ymax = 8;

umax = 10.0; % max accel per axis (hard limit in deterministic MPC) TODO 3.0 causes problem i dont understand

% Chance constraints
beta_state = 0.05;   % probability of violating state bounds (per step, per axis, two-sided)
beta_input = 0.005;   % probability of violating input bounds (per step, per axis, two-sided)

% Chance constraints - 
p_state = 1 - beta_state/2; % constraint is split into two for x and y TODO bunu tam anlamadım
k_state = sqrt(2) * erfinv(2 * p_state - 1);  %inverse CDF

p_input = 1 - beta_input/2;
k_input = sqrt(2) * erfinv(2 * p_input - 1);

% Cost weights
Qpos = 100; % position tracking penalty
Qvel = 2; % velocity penalty
Ru = 0.1; % input penalty
Wx_single = diag([Qpos Qpos Qvel Qvel]); % stage weight on x
Wu_single = Ru*eye(nu); % stage weight on u
WN = 10*Wx_single; % terminal weight

% Waypoints
waypoints = [ 1.0  1.0;
              2.0  6.0;
              6.0  6.0;
              9.5  2.0;
              13.5 6.5;
              15.0 5.0 ];

wp_tol = 0.1; % switch waypoint when within this distance

% Initial true state and initial estimate
x_true0 = [waypoints(1,1); waypoints(1,2); 0; 0];
xhat0   = x_true0 + [0.05; -0.05; 0.1; -0.1];

%% ------------------- Precompute MPC prediction matrices -------------------

% G = [0          0    ...    0
%      B          0    ...    0
%      AB         B    ...    0
%      .          .           .
%      A^(N-1)B   A^(N-2)B    B ]

% H = [I
%      A
%      .
%      A^N]

H = zeros(nx*(N+1), nx);
G = zeros(nx*(N+1), nu*N);

H(1:nx, :) = eye(nx); 

for i = 1:N
    % H matrix : A^i
    H(i*nx+1:(i+1)*nx, :) = A^i;
    
    for j = 1:i
        % G matrix: A^(i-j) * B
        G(i*nx+1:(i+1)*nx, (j-1)*nu+1:j*nu) = (A^(i-j)) * B;
    end
end

% Total state cost matrix
Wx = zeros((N+1)*nx, (N+1)*nx);

% Add wieghts from x_1 to x_{N-1}
for k = 2:N
    idx = (k-1)*nx + (1:nx);
    Wx(idx, idx) = Wx_single;
end

% Add terminal weight x_N
Wx(N*nx+1:end, N*nx+1:end) = WN;

% Total input weights
Wu = kron(eye(N), Wu_single);

%% ------------------- Monte Carlo simulation -------------------

% QP settings
optsQP = optimoptions('quadprog','Display','off');

% J = U'*(G'*Wx*G + Wu)*U + ...
Hessian_qp = G' * Wx * G + Wu; % H in qp
Hessian_qp = (Hessian_qp + Hessian_qp') / 2; % for symmetry

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

        r = waypoints(wp_idx,:)';

        % Check if we should switch to the NEXT centroid
        if wp_idx < size(waypoints,1)
            if norm(x_true(1:2) - r) < wp_tol 
                wp_idx = wp_idx + 1;
                r = waypoints(wp_idx,:)';
            end
        else
            % For the FINAL waypoint, use a tiny tolerance to ensure it stops
            if norm(x_true(1:2) - r) < 0.01 
            end
        end
    
        % Local reference stack
        Xref = repmat([r;0;0], N+1, 1);
        Xref(1:nx) = xhat;    % make x0 reference equal current state

        
        % Future predicted covariance
        Ppred = zeros(nx,nx,N+1);
        Ppred(:,:,1) = P; % x0 covariance
        Ptmp = P;
        for i = 1:N
            Ptmp = A*Ptmp*A' + Q;
            Ppred(:,:,i+1) = Ptmp;
        end
        
        % Chance Constraints
        Aineq = [];
        bineq = [];

        % Set the constraints for the predicted states from i=1 to N (x1...xN).
        for i = 1:N
            % Extract the corresponding px and py indices within the X stack 
            % Starts from x1
            row_px = i*nx + 1;
            row_py = i*nx + 2;
            
            % Matrix G represents the influence of the control inputs
            G_px = G(row_px, :);
            G_py = G(row_py, :);
            
            % Matrix H represents the influence of the initial state x0
            mu_px_aff = H(row_px,:)*xhat;
            mu_py_aff = H(row_py,:)*xhat;
            
            % Tightening 
            sig_px = sqrt(max(Ppred(1,1,i+1),0));
            sig_py = sqrt(max(Ppred(2,2,i+1),0));
            
            lb_px = xmin + k_state*sig_px;
            ub_px = xmax - k_state*sig_px;
            lb_py = ymin + k_state*sig_py;
            ub_py = ymax - k_state*sig_py;
            
            % Aineq*U <= bineq 
            Aineq = [Aineq;  G_px; -G_px;  G_py; -G_py];
            bineq = [bineq;  ub_px - mu_px_aff; -(lb_px - mu_px_aff); ub_py - mu_py_aff; -(lb_py - mu_py_aff)];
        end
        
        % Input constraints
        sig_u = sqrt(Qu(1,1)); % standard deviation of input noise
        umax_tightened = umax - k_input * sig_u;
        
        U_lb = -umax_tightened * ones(nu*N, 1);
        U_ub =  umax_tightened * ones(nu*N, 1);
        
        % Linear Term f: 2 * x0' * H' * Wx * G
        % For reference tracking based on error: 2 * (H*xhat - Xref)' * Wx * G
        f_qp = 2 * (H * xhat - Xref)' * Wx * G;
        f_qp = f_qp'; 
        
        % QP solution
        [Uopt, ~, exitflag] = quadprog(Hessian_qp, f_qp, Aineq, bineq, [], [], U_lb, U_ub, [], optsQP);
        
        if exitflag <= 0 || isempty(Uopt)
            u = zeros(nu,1); % If no solution
        else
            u = Uopt(1:nu);  % Apply first input
        end
        
        % Apply to true system
        w = L_Q * randn(nx, 1);
        x_true = A*x_true + B*u + w;
        
        % Measurement
        v = L_R * randn(ny, 1);
        y = C*x_true + v;
        
        % Kalman filter update
        xpred = A*xhat + B*u;
        PpredKF = A*P*A' + Q;
        S_kf = C*PpredKF*C' + R;
        K = PpredKF*C' / S_kf;
        xhat = xpred + K*(y - C*xpred);
        P = (eye(nx) - K*C)*PpredKF;
        
        % Record
        U_hist(:,k) = u;
        X_hist(:,k) = x_true;
        Xh_hist(:,k)= xhat;
       
    end
    traj_all{mc} = X_hist;
end

%% ------------------- Plot -------------------
figure('Color','w'); 
hold on; 
axis equal; 
grid on;
xlabel('x');
ylabel('y');
title('SMPC Monte Carlo trajectories');

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