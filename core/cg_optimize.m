function result = cg_optimize(prob, D0, S0, method, cgopt)
% CG_OPTIMIZE  Simultaneous recovery of D_i(x,y) and S_i(x,y) by nonlinear CG.
%
%   method : 'FR' (Fletcher-Reeves), 'PR' (Polak-Ribiere+), 'PB' (Powell-Beale),
%            'DY' (Dai-Yuan)
%   D and S are updated together along one search direction with one step
%   length. The variables are scaled (D/Dref, S/Sref) so that both have O(1)
%   size; the CG coefficient beta is computed separately for the D block and
%   the S block (Eqs. (C.1)-(C.5) of the Supplementary Material of the article).
%
%   Line search: 'analytical' (default) = analytical step length from the
%   linearised (sensitivity) problem, as in the classical CGM, followed by
%   backtracking; 'armijo' = heuristic first step + Armijo backtracking.
%   In both cases a step is accepted only if J decreases sufficiently; if no
%   step works the direction is reset to steepest descent once, and the run
%   stops if that also fails.
%   Stopping: (1) projected gradient norm / initial norm < tol_rel,
%             (2) relative decrease of J < ftol for stag_iter consecutive iterations,
%             (3) maxIter, (4) line-search failure. The reason is stored.

    t_start = tic;
    if nargin < 5, cgopt = struct(); end
    cgopt = set_default(cgopt, 'maxIter', 100);
    cgopt = set_default(cgopt, 'tol_rel', 1e-4);
    cgopt = set_default(cgopt, 'ftol', 1e-7);
    cgopt = set_default(cgopt, 'stag_iter', 10);
    cgopt = set_default(cgopt, 'c1', 1e-4);
    cgopt = set_default(cgopt, 'max_bt', 30);
    cgopt = set_default(cgopt, 'step0', 0.05);      % fallback first trial: max change = 5% of max|theta|
    cgopt = set_default(cgopt, 'line_search', 'analytical');  % 'analytical' (+ backtracking) or 'armijo'
    cgopt = set_default(cgopt, 'verbose', true);
    cgopt = set_default(cgopt, 'Cpde', {});          % fields for the PDE-residual (forward error)
    cgopt = set_default(cgopt, 'fixD', false);       % true: keep D fixed, estimate S only (sequential recovery)
    cgopt = set_default(cgopt, 'fixS', false);       % true: keep S fixed, estimate D only

    nsp = prob.nsp;  N = prob.Ny*prob.Nx;  dA = prob.dx*prob.dy;
    nD = nsp*N;
    ip = @(a, b) sum(a.*b)*dA;
    lo = [prob.Dmin/prob.Dref*ones(nD, 1); -inf(nD, 1)];
    hi = [prob.Dmax/prob.Dref*ones(nD, 1);  inf(nD, 1)];
    proj = @(t) min(max(t, lo), hi);
    iD = 1:nD;  iS = nD+1:2*nD;

    th = proj(theta_pack(D0, S0, prob));
    [J, g, info] = eval_theta(th, prob);
    pg = projected_grad(th, g, lo, hi);
    pg = fix_blocks(pg, iD, iS, cgopt);
    gnorm0 = sqrt(ip(pg, pg));

    H = struct('J', [], 'Jdata', [], 'Jreg', [], 'gnorm', [], 'betaD', [], 'betaS', [], ...
               'fwd', [], 'inv', [], 'step', []);
    d_prev = [];  g_prev = [];  s_prev = [];  stag = 0;  restarts = 0;
    stop_reason = 'max_iter';
    if cgopt.verbose
        fprintf('[%s] k=0  J=%.6e (data %.4e, reg %.4e)  |g|=%.3e\n', method, J, info.Jdata, info.Jreg, gnorm0);
    end

    for k = 1:cgopt.maxIter
        gnorm = sqrt(ip(pg, pg));
        H.J(k) = J;  H.Jdata(k) = info.Jdata;  H.Jreg(k) = info.Jreg;  H.gnorm(k) = gnorm;
        H.fwd(k) = forward_error(th, prob, cgopt.Cpde);
        if gnorm/gnorm0 < cgopt.tol_rel
            stop_reason = 'gradient_tolerance';  H.betaD(k) = NaN;  H.betaS(k) = NaN;  H.inv(k) = last_inv(H);  H.step(k) = 0;
            break;
        end

        % ---- search direction ----
        if isempty(d_prev)
            d = -pg;  bD = 0;  bS = 0;
        else
            [bD, rD] = cg_beta(method, pg(iD), g_prev(iD), d_prev(iD), ip);
            [bS, rS] = cg_beta(method, pg(iS), g_prev(iS), d_prev(iS), ip);
            if cgopt.fixD, bD = 0; rD = false; end
            if cgopt.fixS, bS = 0; rS = false; end
            d = -pg;
            if ~rD, d(iD) = -pg(iD) + bD*d_prev(iD); else, bD = 0; restarts = restarts + 1; end
            if ~rS, d(iS) = -pg(iS) + bS*d_prev(iS); else, bS = 0; restarts = restarts + 1; end
            if ~cgopt.fixD && ip(d(iD), pg(iD)) >= 0, d(iD) = -pg(iD); bD = 0; restarts = restarts + 1; end
            if ~cgopt.fixS && ip(d(iS), pg(iS)) >= 0, d(iS) = -pg(iS); bS = 0; restarts = restarts + 1; end
        end
        H.betaD(k) = bD;  H.betaS(k) = bS;

        % ---- initial trial step ----
        if strcmpi(cgopt.line_search, 'analytical')
            s = analytical_step(th, g, d, prob, ip);        % exact minimiser of the linearised J
        elseif isempty(s_prev)
            s = cgopt.step0*max(abs(th))/max(abs(d));
        else
            s = s_prev*ip(g_prev, d_prev)/ip(pg, d);       % Nocedal-Wright heuristic
            s = min(max(s, 0.1*s_prev), 10*s_prev);
        end
        if ~(isfinite(s) && s > 0), s = cgopt.step0*max(abs(th))/max(abs(d)); end

        % ---- Armijo backtracking (projected) ----
        [ok, th_new, J_new, s] = armijo(th, J, g, d, s, prob, proj, cgopt, ip);
        if ~ok && (bD ~= 0 || bS ~= 0)
            d = -pg;  restarts = restarts + 1;
            if strcmpi(cgopt.line_search, 'analytical'), s = analytical_step(th, g, d, prob, ip); else, s = NaN; end
            if ~(isfinite(s) && s > 0), s = cgopt.step0*max(abs(th))/max(abs(d)); end
            [ok, th_new, J_new, s] = armijo(th, J, g, d, s, prob, proj, cgopt, ip);
        end
        if ~ok
            stop_reason = 'line_search_failed';  H.inv(k) = last_inv(H);  H.step(k) = 0;
            break;
        end
        H.step(k) = s;
        H.inv(k) = sqrt(ip(th_new - th, th_new - th))/sqrt(ip(th_new, th_new));

        if (J - J_new)/abs(J) < cgopt.ftol, stag = stag + 1; else, stag = 0; end

        g_prev = pg;  d_prev = d;  s_prev = s;
        th = th_new;
        [J, g, info] = eval_theta(th, prob);
        pg = projected_grad(th, g, lo, hi);
        pg = fix_blocks(pg, iD, iS, cgopt);

        if cgopt.verbose
            fprintf('[%s] k=%d  J=%.6e (data %.4e, reg %.4e)  |g|/|g0|=%.3e  step=%.2e  bD=%.2f bS=%.2f\n', ...
                method, k, J, info.Jdata, info.Jreg, sqrt(ip(pg, pg))/gnorm0, s, bD, bS);
        end
        if stag >= cgopt.stag_iter
            stop_reason = 'stagnation';
            k = k + 1;
            H.J(k) = J;  H.Jdata(k) = info.Jdata;  H.Jreg(k) = info.Jreg;  H.gnorm(k) = sqrt(ip(pg, pg));
            H.fwd(k) = forward_error(th, prob, cgopt.Cpde);  H.betaD(k) = NaN;  H.betaS(k) = NaN;
            H.inv(k) = last_inv(H);  H.step(k) = 0;
            break;
        end
    end
    if strcmp(stop_reason, 'max_iter')
        % record the final state reached after the last update
        k = numel(H.J) + 1;
        H.J(k) = J;  H.Jdata(k) = info.Jdata;  H.Jreg(k) = info.Jreg;  H.gnorm(k) = sqrt(ip(pg, pg));
        H.fwd(k) = forward_error(th, prob, cgopt.Cpde);  H.betaD(k) = NaN;  H.betaS(k) = NaN;
        H.inv(k) = last_inv(H);  H.step(k) = 0;
    end

    [D, S] = theta_unpack(th, prob);
    result = struct();
    result.cg_method = method;
    for s_ = 1:nsp
        result.(sprintf('D%d', s_)) = D{s_};
        result.(sprintf('S%d', s_)) = S{s_};
        result.(sprintf('D%d_phys', s_)) = D{s_};      % already mm^2/s
        result.(sprintf('S%d_phys', s_)) = S{s_};      % already mol/(L s)
    end
    result.C_model = info.C;
    result.J_history = H.J;  result.J_data_history = H.Jdata;  result.J_reg_history = H.Jreg;
    result.grad_norm_history = H.gnorm;
    result.beta_D_history = H.betaD;  result.beta_S_history = H.betaS;
    result.forward_error_history = H.fwd;  result.inverse_error_history = H.inv;
    result.step_history = H.step;
    result.num_iter = numel(H.J);              % recorded points; point 1 = initial guess
    result.num_updates = numel(H.J) - 1;       % accepted parameter updates
    [result.best_cost, ib] = min(H.J);
    result.best_iter = ib;
    result.final_cost = H.J(end);
    result.final_grad_norm = H.gnorm(end);
    result.initial_cost = H.J(1);
    result.cg_restart_count = restarts;
    result.stop_reason = stop_reason;
    result.total_time = toc(t_start);
    if cgopt.verbose
        fprintf('[%s] stop: %s after %d updates, J %.6e -> %.6e (%.2f%% lower), %.1f s\n', ...
            method, stop_reason, result.num_updates, H.J(1), H.J(end), 100*(1 - H.J(end)/H.J(1)), result.total_time);
    end
end

% ------------------------------------------------------------------------
function s = analytical_step(th, g, d, prob, ip)
% Analytical step length of the conjugate gradient method (Orlande-type CGM):
% the step that minimises J along d when the concentration is linearised,
%     C(theta + s d) ~ C + s dC,
% where dC solves the sensitivity problem  A(D) dC = dS - (dA/dD . dD) C
% (zero on the boundary nodes). Then
%     s = - <g, d> / ( sum W dC^2 dxdy + alphaD <d_D, K d_D> + alphaS <d_S, K d_S> ).
    [D, S] = theta_unpack(th, prob);
    [dD, dS] = theta_unpack(d, prob);              % physical directions (scaled back)
    dA = prob.dx*prob.dy;  N = prob.Ny*prob.Nx;  nsp = prob.nsp;
    A0 = adr_assemble(zeros(prob.Ny, prob.Nx), prob);   % advection + boundary rows only
    curv = 0;
    for s_ = 1:nsp
        A  = adr_assemble(D{s_}, prob);
        b  = S{s_}(:);  b(prob.bnd(:)) = prob.gbc{s_};
        C  = A \ b;
        dAC = (adr_assemble(dD{s_}, prob) - A0)*C;     % (dA/dD . dD) C, interior rows only
        rhs = dS{s_}(:) - dAC;  rhs(prob.bnd(:)) = 0;
        dC = A \ rhs;
        curv = curv + sum(prob.W{s_}(:).*dC.^2)*dA;
        tD = d((s_-1)*N + (1:N));  tS = d(nsp*N + (s_-1)*N + (1:N));
        curv = curv + prob.alphaD*(tD.'*(prob.K*tD))*dA + prob.alphaS*(tS.'*(prob.K*tS))*dA;
    end
    s = -ip(g, d)/max(curv, realmin);
end

function [ok, th_new, J_new, s] = armijo(th, J, g, d, s, prob, proj, cgopt, ip)
    ok = false;  th_new = th;  J_new = J;
    for bt = 1:cgopt.max_bt
        th_t = proj(th + s*d);
        J_t = eval_theta(th_t, prob, false);
        if J_t <= J + cgopt.c1*ip(g, th_t - th) && J_t < J
            ok = true;  th_new = th_t;  J_new = J_t;  return;
        end
        s = 0.5*s;
    end
end

function [beta, restart] = cg_beta(method, g, gp, dp, ip)
    restart = false;
    y = g - gp;
    switch upper(method)
        case 'FR'
            beta = ip(g, g)/ip(gp, gp);
        case 'PR'
            beta = max(0, ip(g, y)/ip(gp, gp));
        case 'PB'
            den = ip(dp, y);
            if abs(ip(g, gp)) >= 0.2*ip(g, g) || den <= 0
                beta = 0;  restart = true;
            else
                beta = ip(g, y)/den;
            end
        case 'DY'
            den = ip(dp, y);
            if den <= 0
                beta = 0;  restart = true;
            else
                beta = ip(g, g)/den;
            end
        otherwise
            error('Unknown CG method %s', method);
    end
    if ~isfinite(beta), beta = 0; restart = true; end
end

function [J, g, info] = eval_theta(th, prob, need_grad)
    if nargin < 3, need_grad = true; end
    [D, S] = theta_unpack(th, prob);
    if need_grad
        [J, gD, gS, info] = objective_gradient(D, S, prob, true);
        g = pack_grad(gD, gS, prob);
    else
        J = objective_gradient(D, S, prob, false);
        g = [];  info = [];
    end
end

function v = last_inv(H)
% inverse error of the last accepted update (used for the final record)
    if isempty(H.inv), v = NaN; else, v = H.inv(end); end
end

function pg = fix_blocks(pg, iD, iS, cgopt)
    if cgopt.fixD, pg(iD) = 0; end
    if cgopt.fixS, pg(iS) = 0; end
end

function pg = projected_grad(th, g, lo, hi)
    pg = g;
    pg(th <= lo & g > 0) = 0;
    pg(th >= hi & g < 0) = 0;
end

function e = forward_error(th, prob, Cpde)
% Eqs. (12)-(13): residual of the discrete Eq. (2) evaluated with the
% given concentration fields (e.g. the smooth fitted fields) and the current D, S.
    if isempty(Cpde), e = NaN; return; end
    [D, S] = theta_unpack(th, prob);
    e2 = 0;
    for s = 1:prob.nsp
        A = adr_assemble(D{s}, prob);
        R = A*Cpde{s}(:) - S{s}(:);
        R(prob.bnd(:)) = 0;
        e2 = e2 + sum(R.^2)*prob.dx*prob.dy;
    end
    e = sqrt(e2);
end

function th = theta_pack(D, S, prob)
    th = [];
    for s = 1:prob.nsp, th = [th; D{s}(:)/prob.Dref]; end %#ok<AGROW>
    for s = 1:prob.nsp, th = [th; S{s}(:)/prob.Sref]; end %#ok<AGROW>
end

function g = pack_grad(gD, gS, prob)
    g = [];
    for s = 1:prob.nsp, g = [g; gD{s}(:)*prob.Dref]; end %#ok<AGROW>
    for s = 1:prob.nsp, g = [g; gS{s}(:)*prob.Sref]; end %#ok<AGROW>
end

function [D, S] = theta_unpack(th, prob)
    N = prob.Ny*prob.Nx;  D = cell(1, prob.nsp);  S = cell(1, prob.nsp);
    for s = 1:prob.nsp
        D{s} = reshape(th((s-1)*N + (1:N)), prob.Ny, prob.Nx)*prob.Dref;
        S{s} = reshape(th(prob.nsp*N + (s-1)*N + (1:N)), prob.Ny, prob.Nx)*prob.Sref;
    end
end

function s = set_default(s, f, v)
    if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end
