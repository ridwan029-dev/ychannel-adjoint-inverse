function LC2 = lcurve_sweep2(x, y, v, C_obs, W, C_bc, opt, D0, S0, alphasD, alphasS, method, cgopt, ...
                             alpha_start, n_cycles, cache, fname)
% LCURVE_SWEEP2  Separate regularization weights alphaD and alphaS by alternating L-curves.
%   The two-parameter problem is reduced to a sequence of one-parameter L-curves:
%     (a) alphaS fixed: sweep alphaD over ALPHASD and take the corner of
%         (data misfit, H1 seminorm of the scaled D fields)  -> new alphaD;
%     (b) alphaD fixed: sweep alphaS over ALPHASS and take the corner of
%         (data misfit, H1 seminorm of the scaled S fields)  -> new alphaS.
%   The cycle starts from alphaD = alphaS = ALPHA_START (e.g. the corner of the joint
%   L-curve) and is repeated at most N_CYCLES times; it stops when neither value changes.
%   Every (alphaD, alphaS) pair is run only once (CACHE = rows [aD aS misfit etaD etaS],
%   e.g. LC.runs from lcurve_sweep, may be passed in).
%   Output: LC2.alphaD, LC2.alphaS, LC2.runs (all runs), LC2.sweeps (one entry per sweep).
    if nargin < 15 || isempty(n_cycles), n_cycles = 2; end
    if nargin < 16, cache = zeros(0, 5); end
    if nargin < 17 || isempty(fname), fname = 'L_curve_separate.png'; end
    cgopt.verbose = false;
    ctx = struct('x', x, 'y', y, 'v', v, 'opt', opt, 'method', method, 'cgopt', cgopt);
    ctx.C_obs = C_obs;  ctx.W = W;  ctx.C_bc = C_bc;  ctx.D0 = D0;  ctx.S0 = S0;
    aD = alpha_start;  aS = alpha_start;
    runs = cache;  sweeps = {};
    for c = 1:n_cycles
        aD_old = aD;  aS_old = aS;
        % (a) sweep alphaD with alphaS fixed
        T = zeros(numel(alphasD), 3);
        for k = 1:numel(alphasD)
            [row, runs] = get_run(ctx, alphasD(k), aS, runs);
            T(k, :) = [alphasD(k), row(3), row(4)];
        end
        aD = pick_corner(T, 'alphaD');
        sweeps{end + 1} = struct('cycle', c, 'which', 'D', 'fixed', aS, 'table', T, 'selected', aD); %#ok<AGROW>
        fprintf('   cycle %d: alphaS = %.1e fixed -> alphaD = %.1e\n', c, aS, aD);
        % (b) sweep alphaS with alphaD fixed
        T = zeros(numel(alphasS), 3);
        for k = 1:numel(alphasS)
            [row, runs] = get_run(ctx, aD, alphasS(k), runs);
            T(k, :) = [alphasS(k), row(3), row(5)];
        end
        aS = pick_corner(T, 'alphaS');
        sweeps{end + 1} = struct('cycle', c, 'which', 'S', 'fixed', aD, 'table', T, 'selected', aS); %#ok<AGROW>
        fprintf('   cycle %d: alphaD = %.1e fixed -> alphaS = %.1e\n', c, aD, aS);
        if aD == aD_old && aS == aS_old, break; end
    end
    LC2 = struct('alphaD', aD, 'alphaS', aS, 'runs', runs, 'cycles', c);
    LC2.sweeps = sweeps;
    fprintf('   -> separate weights: alphaD = %.1e, alphaS = %.1e (%d cycle(s), %d inversions in total)\n', ...
        aD, aS, c, size(runs, 1));
    plot_separate(sweeps, aD, aS, fname);
end


function [row, runs] = get_run(ctx, a_d, a_s, runs)
% one inversion for (alphaD, alphaS), or the cached result if this pair was run already
    if ~isempty(runs)
        i = find(abs(log10(runs(:, 1)) - log10(a_d)) < 1e-9 & abs(log10(runs(:, 2)) - log10(a_s)) < 1e-9, 1);
        if ~isempty(i), row = runs(i, :); return; end
    end
    o = ctx.opt;  o.alphaD = a_d;  o.alphaS = a_s;
    prob = make_problem(ctx.x, ctx.y, ctx.v, ctx.C_obs, ctx.W, ctx.C_bc, o);
    r = cg_optimize(prob, ctx.D0, ctx.S0, ctx.method, ctx.cgopt);
    [eD, eS] = scaled_seminorms(r, prob);
    row = [a_d, a_s, r.J_data_history(end), eD, eS];
    runs(end + 1, :) = row;
    fprintf('   alphaD = %.1e, alphaS = %.1e: data misfit %.4e, seminorm D %.4e, S %.4e (%s)\n', ...
        a_d, a_s, row(3), eD, eS, r.stop_reason);
end

function a = pick_corner(T, label)
    [a, ~, at_edge] = lcurve_choice(T, 1);
    if isnan(a), [~, i] = min(T(:, 2)); a = T(i, 1); end
    if at_edge
        warning(['The L-curve of %s bends more strongly at an end of the tested range than at the ', ...
                 'selected corner (%.1e); extend the range to confirm the corner.'], label, a);
    end
end

function plot_separate(sweeps, aD, aS, fname)
    try
        figH = figure('Color', 'w', 'Name', 'Separate L-curves', 'Position', [60 60 1100 420]);
        last = {[], []};
        for k = 1:numel(sweeps)
            if sweeps{k}.which == 'D', last{1} = sweeps{k}; else, last{2} = sweeps{k}; end
        end
        lab = {'H^1 seminorm of scaled D', 'H^1 seminorm of scaled S'};
        sel = [aD, aS];  nm = {'\alpha_D', '\alpha_S'};  fx = {'\alpha_S', '\alpha_D'};
        for p = 1:2
            T = last{p}.table;
            subplot(1, 2, p);
            loglog(T(:, 2), T(:, 3), 'ko-', 'LineWidth', 1.5); hold on;
            i = find(T(:, 1) == sel(p), 1);
            loglog(T(i, 2), T(i, 3), 'rs', 'MarkerSize', 12, 'LineWidth', 2); hold off; grid on;
            text(T(:, 2), T(:, 3), arrayfun(@(a) sprintf('  %.0e', a), T(:, 1), 'UniformOutput', false));
            xlabel('data misfit'); ylabel(lab{p});
            title(sprintf('%s sweep (%s = %.0e): selected %.0e', nm{p}, fx{p}, last{p}.fixed, sel(p)));
        end
        lc_save_png(figH, fname);
    catch ME
        fprintf('   (separate L-curve figure not saved: %s)\n', ME.message);
    end
end
