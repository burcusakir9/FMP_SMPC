function out = simulateRun(chain, cfg, run)
%SIMULATERUN  One vessel run through an RSC funnel chain.
%
%   out = simulateRun(chain, cfg, run)
%
% run fields (defaults in brackets):
%   law     'durmaz' | 'ms' | 'msR'      speed law (see thesisConfig.m)
%   U_m     mission speed [m/s]          (ignored by 'durmaz')
%   filter  ['none'] | 'cbf' | 'hocbf'   safety filter; 'rcbf' | 'rhocbf':
%           robust versions for a current up to cfg.Vb
%   plant   ['dynamic'] | 'unicycle'     vessel model or ideal kinematic unicycle
%   dist    struct of disturbances, all unknown to controller and filters:
%             Vc [0;0]        constant current [m/s], world frame
%             ins_pos [0]     INS position noise std [m]
%             ins_psi [0]     INS heading noise std [rad]
%             ins_vel [0]     INS noise std on u, v, r
%             act_gain [1;1]  applied/commanded thrust, [left; right]
%             act_noise [0]   thrust noise std [N]
%             seed [1]        noise seed
%   keepLog [false]           return the time histories (0.1 s)
%
% Control loop (every dt):
%   1. measured state xm = true state + INS noise
%   2. active funnel = highest-index path funnel containing xm (if none,
%      the previous one is kept)
%   3. Durmaz2024 law with speed s(rho): u = s*cos(alpha),
%      w = Ka*alpha + (u/rho)*sin(alpha)
%   4. 'cbf'  : kinematic CBF (u) / HOCBF (w) filter on (u, w)
%   5. PI loops on u and r (conditional-integration anti-windup) and
%      thrust allocation F_L,R = X/2 +- N/(2d), saturated to the limits
%   6. 'hocbf': second-order HOCBF filter on the thrust (vessel model)
%   7. plant: thrust efficiency/noise, RK4 of the vessel model + current
%
% Metrics (true state):
%   reached, stuck, T           goal reached / no progress for cfg.Tstuck s
%   exitAct                     max distance outside the active funnel [m]; the
%                               active funnel for the metrics is chosen with the same
%                               rule from the TRUE position (with INS noise the
%                               controller's funnel switches early on the measured one)
%   exitChain, tOutChain        max distance / time outside every path funnel
%   collision                   hull reference point inside an obstacle
%   uMean, uCV                  surge speed mean and std/mean (intermediate funnels)
%   filtPct, infeasPct          steps the filter changed the command /
%                               steps its condition was unsatisfiable [%]
%   funnelExit (1 x n)          max exit while funnel k was active (NaN: never)

    d0 = struct('filter', 'none', 'plant', 'dynamic', 'dist', struct(), 'keepLog', false, 'U_m', NaN);
    for f = fieldnames(d0)', if ~isfield(run, f{1}), run.(f{1}) = d0.(f{1}); end, end
    D = struct('Vc', [0; 0], 'ins_pos', 0, 'ins_psi', 0, 'ins_vel', 0, ...
               'act_gain', [1; 1], 'act_noise', 0, 'seed', 1);
    for f = fieldnames(run.dist)', D.(f{1}) = run.dist.(f{1}); end

    C = chain.C; R = chain.R; n = numel(R);
    vs = cfg.vessel; model = cfg.model; dt = cfg.dt;
    dyn = strcmp(run.plant, 'dynamic');

    % Simulation horizon from the path length and the nominal speed
    Lpath = sum(vecnorm(diff([chain.q_start, C], 1, 2)));
    if strcmp(run.law, 'durmaz'), Unom = cfg.Kv*mean(R); else, Unom = run.U_m; end
    Tmax = max(500, cfg.Tfactor*Lpath/Unom);
    N = ceil(Tmax/dt);

    KF = cfg.KF; KF.u_lim = max(1, run.U_m*~strcmp(run.law, 'durmaz')); KF.qp = cfg.qp;
    HF = cfg.HF; HF.qp = cfg.qp;
    robust = any(strcmp(run.filter, {'rcbf', 'rhocbf'}));
    KF.Vb = robust*cfg.Vb; HF.Vb = robust*cfg.Vb;
    useKF = any(strcmp(run.filter, {'cbf', 'rcbf'}));
    useHF = any(strcmp(run.filter, {'hocbf', 'rhocbf'}));
    PI = cfg.PI;
    Flim = [vs.Fmin, vs.Fmax];

    rng(D.seed);
    g = C(:, min(2, n)) - chain.q_start;
    x = [chain.q_start; atan2(g(2), g(1)); 0; 0; 0; 0; 0];   % at rest, facing the 2nd funnel
    int_u = 0; int_r = 0;
    kAct = 1; kTrue = 1; kMax = 1; tLast = 0;

    nLog = floor(N/cfg.logEvery) + 1; iLog = 0;
    L = struct('t', nan(1,nLog), 'X', nan(1,nLog), 'Y', nan(1,nLog), 'u', nan(1,nLog), ...
               'r', nan(1,nLog), 'h', nan(1,nLog), 'k', nan(1,nLog), 'filt', false(1,nLog));
    funnelExit = nan(1, n);
    exitAct = 0; exitChain = 0; tOut = 0;
    nFilt = 0; nInfeas = 0; nCtrl = 0;
    uSum = 0; u2Sum = 0; uN = 0;
    reached = false; stuck = false; t = 0;

    for step = 1:N
        t = (step - 1)*dt;

        % 1. Measured state
        xm = x;
        xm(1:2) = x(1:2) + D.ins_pos*randn(2,1);
        xm(3)   = wrapPi(x(3) + D.ins_psi*randn);
        xm(4:6) = x(4:6) + D.ins_vel*randn(3,1);

        % 2. Active funnel (measured position)
        inside = find(vecnorm(C - xm(1:2)) < R, 1, 'last');
        if ~isempty(inside), kAct = inside; end
        if kAct > kMax, kMax = kAct; tLast = t; end
        ctr = C(:, kAct); Rk = R(kAct); isGoal = (kAct == n);

        % Metrics on the true state
        inT = find(vecnorm(C - x(1:2)) < R, 1, 'last');
        if ~isempty(inT), kTrue = inT; end
        e1 = max(0, norm(x(1:2) - C(:, kTrue)) - R(kTrue));
        e2 = max(0, min(vecnorm(C - x(1:2)) - R));
        exitAct = max(exitAct, e1); exitChain = max(exitChain, e2);
        tOut = tOut + dt*(e2 > 0);
        funnelExit(kTrue) = max([funnelExit(kTrue), e1], [], 'omitnan');
        if kAct < n, uSum = uSum + x(4); u2Sum = u2Sum + x(4)^2; uN = uN + 1; end

        if mod(step - 1, cfg.logEvery) == 0
            iLog = iLog + 1;
            L.t(iLog) = t; L.X(iLog) = x(1); L.Y(iLog) = x(2); L.u(iLog) = x(4);
            L.r(iLog) = x(6); L.h(iLog) = R(kTrue) - norm(x(1:2) - C(:, kTrue)); L.k(iLog) = kTrue;
        end

        if norm(x(1:2) - chain.q_goal) <= cfg.goal_tol, reached = true; break; end
        if t - tLast > cfg.Tstuck, stuck = true; break; end

        % 3. Durmaz2024 law with the speed law s(rho)
        [rho, alpha] = polar(xm, ctr);
        switch run.law
            case 'durmaz', U = NaN;
            case 'ms',     U = run.U_m;
            case 'msR',    U = min(run.U_m, Rk/cfg.T_R);
        end
        if strcmp(run.law, 'durmaz')
            s = 2*cfg.Kv*rho;
        elseif isGoal
            s = max(U*tanh(2*cfg.Kv*rho/U), min(U, cfg.U_goalMin));
        else
            s = U;
        end
        if rho <= cfg.rho_tol
            un = 0; wn = 0;
        else
            un = s*cos(alpha);
            wn = cfg.Ka*alpha + (un/rho)*sin(alpha);
        end
        u_ref = un; w_ref = wn; changed = false;

        % 4. Kinematic CBF/HOCBF filter on the references
        if useKF && rho > cfg.rho_tol
            [u_ref, w_ref, infeas] = kinFilter(xm, un, wn, ctr, Rk, KF);
            changed = abs(u_ref - un) > 1e-4 || abs(w_ref - wn) > 1e-4;
            nInfeas = nInfeas + infeas; nCtrl = nCtrl + 1;
        end

        if ~dyn
            % Ideal kinematic unicycle: the references are realized exactly
            x(1:3) = x(1:3) + dt*[u_ref*cos(x(3)) + D.Vc(1); u_ref*sin(x(3)) + D.Vc(2); w_ref];
            x(3) = wrapPi(x(3)); x(4) = u_ref; x(5) = 0; x(6) = w_ref;
        else
            % 5. PI loops (conditional integration) + allocation
            e_u = u_ref - xm(4); e_r = w_ref - xm(6);
            iu = min(max(int_u + e_u*dt, -PI.Ispeed_max), PI.Ispeed_max);
            ir = min(max(int_r + e_r*dt, -PI.Iyaw_max), PI.Iyaw_max);
            tauX = PI.P_speed*e_u + PI.I_speed*iu;
            tauN = PI.P_yaw*e_r + PI.I_yaw*ir;
            FLc = tauX/2 + tauN/(2*vs.d); FRc = tauX/2 - tauN/(2*vs.d);
            F = [min(max(FLc, Flim(1)), Flim(2)); min(max(FRc, Flim(1)), Flim(2))];
            if F(1) == FLc && F(2) == FRc, int_u = iu; int_r = ir; end

            % 6. Dynamic HOCBF filter on the thrust
            if useHF && rho > cfg.rho_tol
                Fn = F;
                [F, infeas] = hocbfFilter(xm, Fn, ctr, Rk, HF, Flim, model, vs);
                changed = norm(F - Fn) > 1e-3;
                nInfeas = nInfeas + infeas; nCtrl = nCtrl + 1;
            end

            % 7. Actuator disturbance and plant (RK4) with the current
            Fa = D.act_gain(:).*F + D.act_noise*randn(2,1);
            cur = [D.Vc(:); zeros(6,1)];
            k1 = model(x, Fa) + cur;
            k2 = model(x + 0.5*dt*k1, Fa) + cur;
            k3 = model(x + 0.5*dt*k2, Fa) + cur;
            k4 = model(x + dt*k3, Fa) + cur;
            x = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
            x(3) = wrapPi(x(3));
        end
        nFilt = nFilt + changed;
        if mod(step - 1, cfg.logEvery) == 0, L.filt(iLog) = changed; end
    end

    % Trim the logs; obstacle check on the logged positions
    for f = fieldnames(L)', L.(f{1}) = L.(f{1})(1:iLog); end
    collision = any(isinterior(chain.obsAll, L.X(:), L.Y(:)));

    out.reached = reached; out.stuck = stuck; out.T = t;
    out.exitAct = exitAct; out.exitChain = exitChain; out.tOutChain = tOut;
    out.collision = collision;
    out.uMean = uSum/max(uN, 1);
    out.uCV = sqrt(max(u2Sum/max(uN, 1) - out.uMean^2, 0))/max(abs(out.uMean), 1e-9);
    out.filtPct = 100*nFilt/step;
    out.infeasPct = 100*nInfeas/max(nCtrl, 1);
    out.funnelExit = funnelExit;
    if run.keepLog, out.log = L; end
end

%% ---------------- Safety filters ----------------

function [u_ref, w_ref, infeas] = kinFilter(x, u_nom, w_nom, ctr, R, KF)
% Kinematic CBF/HOCBF filter on nu = [u; w] (compareFunnelChain.m, CBF.m).
% Model xi_dot = R(psi)*[0; v] + g(xi)*nu with the measured sway v as drift;
% barrier b = R - rho:
%   (1) CBF on u   : Lf b + Lg b*nu + k1*b >= 0
%   (2) HOCBF on w : Lf2 b + LgLf b*nu + (k1+k2)*b_dot + k1*k2*b >= 0
%       (surge frozen at its measured value, u_dot and v_dot neglected)
% QP: min ||nu - nu_nom||^2 + Wd*(d1^2 + d2^2) with slacks d >= 0.
    [rho, alpha] = polar(x, ctr);
    u = x(4); v = x(5);
    wt = -u*sin(alpha) + v*cos(alpha);            % tangential velocity
    b = R - rho;
    Lfb = v*sin(alpha); Lgb = [cos(alpha), 0];
    bdot = Lfb + Lgb(1)*u;
    Lf2b = -wt^2/rho; LgLfb = [0, -wt];

    A = [-Lgb; -LgLfb];
    % Robust (KF.Vb > 0): b_dot is replaced by its worst case b_dot - Vb
    c = [Lfb - KF.Vb + KF.k1*b; Lf2b + (KF.k1 + KF.k2)*(bdot - KF.Vb) + KF.k1*KF.k2*b];
    lb = [-KF.u_lim; -KF.w_lim]; ub = [KF.u_lim; KF.w_lim];
    nu_nom = [u_nom; w_nom];

    best = sum(min(A.*lb', A.*ub'), 2);           % least A*nu over the box, per row
    infeas = any(best > c);
    u_ref = u_nom; w_ref = w_nom;
    if all(A*nu_nom <= c), return; end

    H = diag([2, 2, 2*KF.slack_w, 2*KF.slack_w]);
    f = [-2*nu_nom; 0; 0];
    [z, ~, flag] = quadprog(H, f, [A, -eye(2)], c, [], [], [lb; 0; 0], [ub; inf; inf], ...
        [min(max(nu_nom, lb), ub); 0; 0], KF.qp);
    if flag > 0, u_ref = z(1); w_ref = z(2); end
end

function [F, infeas] = hocbfFilter(x, F_nom, ctr, R, HF, Flim, model, vs)
% Second-order HOCBF on the thrust F = [F_L; F_R] (compareFunnelChain.m, HOCBF.m).
% b = R - rho has relative degree 2 in F:
%   b_ddot = Lf2b + LgLfb*F,  Lf2b = cos(a)*d1 + sin(a)*d2 - wt^2/rho - wt*r,
%   LgLfb = cos(a)*G(1,:) + sin(a)*G(2,:),  d = M^-1(-C nu - D nu),  G = M^-1 B
% (d from the vessel model with zero thrust).
% QP: min ||F - F_nom||^2 + Wd*delta^2
%     s.t. b_ddot + (k1+k2)*b_dot + k1*k2*b >= -delta, delta >= 0, Flim.
    [rho, alpha] = polar(x, ctr);
    nu = x(4:6); u = nu(1); v = nu(2); r = nu(3);
    wt = -u*sin(alpha) + v*cos(alpha);
    xd = model([0; 0; 0; nu; 0; 0], [0; 0]); dacc = xd(4:6);
    G = vs.M \ vs.B_prop;

    b = R - rho;
    Lfb = u*cos(alpha) + v*sin(alpha);
    Lf2b = cos(alpha)*dacc(1) + sin(alpha)*dacc(2) - wt^2/rho - wt*r;
    LgLfb = cos(alpha)*G(1,:) + sin(alpha)*G(2,:);

    A = -LgLfb;
    c = Lf2b + (HF.k1 + HF.k2)*(Lfb - HF.Vb) + HF.k1*HF.k2*b;   % robust if HF.Vb > 0
    infeas = sum(min(A*Flim(1), A*Flim(2))) > c;
    F = F_nom;
    if A*F_nom <= c, return; end

    H = diag([2, 2, 2*HF.slack_w]);
    f = [-2*F_nom; 0];
    [z, ~, flag] = quadprog(H, f, [A, -1], c, [], [], [Flim(1); Flim(1); 0], ...
        [Flim(2); Flim(2); inf], [F_nom; 0], HF.qp);
    if flag > 0, F = z(1:2); end
end

%% ---------------- Helpers ----------------

function [rho, alpha] = polar(x, ctr)
    ex = x(1) - ctr(1); ey = x(2) - ctr(2);
    rho = hypot(ex, ey);
    alpha = wrapPi(atan2(-ey, -ex) - x(3));
end

function a = wrapPi(a)
    a = mod(a + pi, 2*pi) - pi;
end
