function [Cfit, P, stats] = fit_profiles(Cobs, x, y, model, opts)
% FIT_PROFILES  Least-squares fit of analytical profiles to the measured maps.
%
%   Cobs  : Ny x Nx measured concentration (mol/L), NaN allowed
%   x     : 1 x Nx transverse coordinate (any unit; same unit for mu and sigma)
%   y     : Ny x 1 streamwise coordinate
%   model : 'erfc_left'  C = (a/2) erfc((x - mu)/s)     (H2SO4, high at x < mu)
%           'erfc_right' C = (a/2) erfc(-(x - mu)/s)    (NaOH,  high at x > mu)
%           'gauss'      C = a exp(-(x - mu)^2/(2 s^2)) (Na2SO4)
%   opts.ydeg  : polynomial degree of the parameters in y (default 1 = linear)
%   opts.rows  : rows used for the row-by-row fits (default: all)
%
% Step 1: fit (a, mu, s) separately for every row y_i (fminsearch, no toolbox).
% Step 2: fit a polynomial in y to each parameter (weighted by the row R^2),
%         and rebuild a smooth field Cfit(x, y).
% stats: RMSE and R^2 of Cfit against Cobs, and the median row-wise R^2.

    if nargin < 5, opts = struct(); end
    if ~isfield(opts, 'ydeg'), opts.ydeg = 1; end
    x = x(:).';  y = y(:);
    [Ny, Nx] = size(Cobs);
    if ~isfield(opts, 'rows'), opts.rows = 1:Ny; end
    f = profile_fun(model);
    span = max(x) - min(x);

    rows = opts.rows(:).';
    prm = nan(numel(rows), 3);  r2row = nan(numel(rows), 1);
    dxs = abs(x(2) - x(1));  smin = dxs/2;  smax = span;
    gmax = max(Cobs(:));  grange = gmax - min(Cobs(:));
    amax = 2*gmax + eps;  xmid = 0.5*(max(x) + min(x));  mur = 0.75*span;
    % bounded parameterisation: a in (0, amax), mu in xmid +/- mur, s in (smin, smax)
    tr  = @(q) [amax./(1 + exp(-q(1))), xmid + mur*tanh(q(2)), smin + (smax - smin)./(1 + exp(-q(3)))];
    itr = @(p) [logit_(p(1)/amax), atanh_(min(max((p(2) - xmid)/mur, -0.999), 0.999)), logit_((p(3) - smin)/(smax - smin))];
    fo = optimset('Display', 'off', 'MaxFunEvals', 4000, 'MaxIter', 4000, 'TolX', 1e-8, 'TolFun', 1e-10);
    prev = [];
    for k = 1:numel(rows)
        c = Cobs(rows(k), :);  ok = ~isnan(c);
        if sum(ok) < 6, continue; end
        if max(c(ok)) - min(c(ok)) < 0.05*grange, continue; end     % too little signal in this row
        xc = x(ok);  cc = c(ok);
        p0 = initial_guess(model, xc, cc, span);
        if ~isempty(prev), p0 = prev; end                    % warm start from the previous row
        obj = @(q) sum((f(tr(q), xc) - cc).^2);
        q = fminsearch(obj, itr(p0), fo);
        % restart from the data-based guess if the warm start ended badly
        if ~isempty(prev)
            q1 = fminsearch(obj, itr(initial_guess(model, xc, cc, span)), fo);
            if obj(q1) < obj(q), q = q1; end
        end
        p = tr(q);
        prm(k, :) = p;  prev = p;
        r2row(k) = 1 - obj(q)/sum((cc - mean(cc)).^2);
    end

    good = all(isfinite(prm), 2);
    w = max(r2row(good), 0) + 1e-3;
    yy = y(rows(good));
    P.coef = zeros(3, opts.ydeg + 1);
    for j = 1:3
        P.coef(j, :) = wpolyfit(yy, prm(good, j), opts.ydeg, w);
    end
    a  = polyval(P.coef(1, :), y);
    mu = polyval(P.coef(2, :), y);
    s  = min(max(polyval(P.coef(3, :), y), smin), smax);
    a  = min(max(a, 0), amax);
    Cfit = zeros(Ny, Nx);
    for i = 1:Ny, Cfit(i, :) = f([a(i), mu(i), s(i)], x); end

    P.a = a;  P.mu = mu;  P.s = s;  P.model = model;
    P.row_params = prm;  P.row_R2 = r2row;  P.rows = rows;
    m = ~isnan(Cobs);
    res = Cfit(m) - Cobs(m);
    stats.RMSE = sqrt(mean(res.^2));
    stats.R2 = 1 - sum(res.^2)/sum((Cobs(m) - mean(Cobs(m))).^2);
    stats.median_row_R2 = median(r2row(good));
end

function v = logit_(t)
    t = min(max(t, 1e-6), 1 - 1e-6);
    v = log(t/(1 - t));
end

function v = atanh_(t)
    v = 0.5*log((1 + t)/(1 - t));
end

function f = profile_fun(model)
    switch lower(model)
        case 'erfc_left',  f = @(p, x) 0.5*p(1)*erfc((x - p(2))/p(3));
        case 'erfc_right', f = @(p, x) 0.5*p(1)*erfc(-(x - p(2))/p(3));
        case 'gauss',      f = @(p, x) p(1)*exp(-(x - p(2)).^2/(2*p(3)^2));
        otherwise, error('unknown model %s', model);
    end
end

function p = initial_guess(model, x, c, span)
    switch lower(model)
        case {'erfc_left', 'erfc_right'}
            a = max(prctile_simple(c, 95), 1e-3);
            [~, i] = min(abs(c - a/2));
            p = [a, x(i), span/10];
        case 'gauss'
            [a, i] = max(c);
            p = [max(a, 1e-3), x(i), span/10];
    end
end

function v = prctile_simple(c, q)
    c = sort(c(:));
    v = c(max(1, round(q/100*numel(c))));
end

function p = wpolyfit(x, y, n, w)
    V = zeros(numel(x), n + 1);
    for k = 0:n, V(:, n + 1 - k) = x(:).^k; end
    sw = sqrt(w(:));
    p = ((V.*sw) \ (y(:).*sw)).';
end
