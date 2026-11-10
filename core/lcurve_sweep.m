function LC = lcurve_sweep(x, y, v, C_obs, W, C_bc, opt, D0, S0, alphas, method, cgopt, noise_sd, rule, fname)
% LCURVE_SWEEP  L-curve for alphaD = alphaS = alpha and the choice of alpha.
%   For every alpha in ALPHAS one inversion is run (METHOD, CGOPT) and the
%   data misfit and the H1 seminorm of the scaled fields are recorded.
%   rule = 'corner'      : alpha at the maximum curvature of the L-curve (default)
%          'discrepancy' : alpha whose data misfit is closest to the value expected
%                          from the measurement noise NOISE_SD (mol/L)
%   LC.alpha is the selected value; LC.table = [alpha, misfit, seminorm].
%   LC.runs = [alphaD, alphaS, misfit, seminorm of D, seminorm of S] (used by lcurve_sweep2).
%   A figure is saved to FNAME (default 'L_curve.png').
    if nargin < 14 || isempty(rule), rule = 'corner'; end
    if nargin < 15 || isempty(fname), fname = 'L_curve.png'; end
    cgopt.verbose = false;
    Lc = zeros(numel(alphas), 3);  runs = zeros(numel(alphas), 5);
    for k = 1:numel(alphas)
        opt.alphaD = alphas(k);  opt.alphaS = alphas(k);
        prob = make_problem(x, y, v, C_obs, W, C_bc, opt);
        r = cg_optimize(prob, D0, S0, method, cgopt);
        Lc(k, :) = [alphas(k), r.J_data_history(end), r.J_reg_history(end)/alphas(k)];
        [eD, eS] = scaled_seminorms(r, prob);
        runs(k, :) = [alphas(k), alphas(k), r.J_data_history(end), eD, eS];
        fprintf('   L-curve alpha = %.1e: data misfit %.4e, seminorm %.4e (%s)\n', ...
            Lc(k, 1), Lc(k, 2), Lc(k, 3), r.stop_reason);
    end
    nvalid = sum(cellfun(@(w) sum(w(:)), W));
    target = 0.5*nvalid*noise_sd^2*(x(2) - x(1))*(y(2) - y(1));
    [a_corner, a_disc, at_edge] = lcurve_choice(Lc, target);
    switch lower(rule)
        case 'discrepancy', a_sel = a_disc;
        otherwise,          a_sel = a_corner;
    end
    if isnan(a_sel), a_sel = a_disc; end
    if strcmpi(rule, 'corner') && at_edge
        warning(['The L-curve bends more strongly at an end of the tested range than at the selected ', ...
                 'corner (%.1e); extend ALPHAS to confirm the corner.'], a_sel);
    elseif ~strcmpi(rule, 'corner') && any(a_sel == [alphas(1), alphas(end)])
        warning('The selected alpha (%.1e) is at the end of the tested range; extend ALPHAS.', a_sel);
    end
    LC = struct('table', Lc, 'runs', runs, 'alpha', a_sel, 'alpha_corner', a_corner, ...
                'alpha_discrepancy', a_disc, 'target_misfit', target, 'rule', rule);
    fprintf('   -> L-curve corner alpha = %.1e, discrepancy alpha = %.1e; using %.1e (%s)\n', ...
        a_corner, a_disc, a_sel, rule);
    try
        figH = figure('Color', 'w', 'Name', 'L-curve');
        loglog(Lc(:, 2), Lc(:, 3), 'ko-', 'LineWidth', 1.5); hold on;
        i = find(Lc(:, 1) == a_sel, 1);
        loglog(Lc(i, 2), Lc(i, 3), 'rs', 'MarkerSize', 12, 'LineWidth', 2); hold off; grid on;
        text(Lc(:, 2), Lc(:, 3), arrayfun(@(a) sprintf('  %.0e', a), Lc(:, 1), 'UniformOutput', false));
        xlabel('data misfit'); ylabel('H^1 seminorm of scaled D and S');
        title(sprintf('L-curve (selected \\alpha = %.0e, %s)', a_sel, rule));
        lc_save_png(figH, fname);
    catch ME
        fprintf('   (L-curve figure not saved: %s)\n', ME.message);
    end
end
