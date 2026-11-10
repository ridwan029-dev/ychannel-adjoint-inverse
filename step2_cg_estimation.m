%STEP2_CG_ESTIMATION  Simultaneous recovery of D_i(x,y) and S_i(x,y) with
%   adjoint-based nonlinear conjugate gradient methods.
%
%   Step 2 of 3. Starting from the estimates of step 1 for the same data set,
%   the script
%     A. fits the profiles, builds the initial guess and sets the
%        regularisation weights (fixed values or L-curves),
%     B. runs the four CG methods (FR, PR, PB, DY) and the sequential
%        recovery (D fixed),
%     C. prints the comparison tables and plots the estimated fields,
%     D. runs the sensitivity analyses and the Monte Carlo uncertainty analysis.
%
%   Method (Sections 3-4 and Appendices A-C)
%    - Forward model, Eq. (2): v dC/dy = d/dx( D dC/dx ) + S, conservative finite
%      volumes in x; solved by marching along the channel (y) from the inlet, with
%      all points of each cross-section solved simultaneously (implicit scheme);
%      boundary values from the fitted fields.
%    - Exact discrete adjoint for the gradients (checked against finite differences).
%    - D and S are updated together along one search direction (scaled variables);
%      beta is computed separately for the D and S blocks (Eqs. (C.1)-(C.5)).
%    - Step length from the linearised (sensitivity) problem with backtracking;
%      a step is accepted only if J decreases.
%    - H1 regularisation; weights fixed or from L-curves (alpha_mode).
%    - Physical units throughout (mm, s, mol/L, mol/(L s)).
%
%   Input    results/<condition>/step1_initial_estimates.mat (from step 1)
%   Output   in results/<condition>/:
%            results_summary.txt    all numbers of the run (tables, sensitivity, Monte Carlo)
%            step2_cg_results.mat   estimated fields of every method (input of step 3)
%            *.png                  figures
%            With lcurve_only = true only step2_lcurve_results.txt/.mat are written.
%   Results  Figs. 3 and 5-11, Tables 1-2, Sections 6.2-6.5 (condition T28);
%            Table 4, Section 6.6 (one run per condition)
%   Usage    Run after step 1 with the same condition. All settings are in the
%            "Settings" block below; the defaults reproduce the article.
%   Requires MATLAB R2020a or newer (base MATLAB only). All functions are local
%            functions at the end of this file; the same functions are also in core/.
%   Runtime  Several hours with the defaults (300 iterations per CG method and
%            the Monte Carlo analysis with N_MC = 20); set run_analyses = false
%            for a shorter run.
%
%   Code accompanying the article
%     R. Abdurrahman, K. Matsunaga, Y. Abe, N. Kakuta, "Simultaneous Regularized
%     Field Recovery of Diffusion Coefficients and Source Terms in Microfluidic
%     Acid-Base Reactions",
%     Chemical Engineering Science (submitted; DOI to be added on publication).
%   Equation, figure and table numbers in the comments refer to this article
%   (equations C.1-C.5: its Supplementary Material).
%
%   Repository : https://github.com/ridwan029-dev/ychannel-adjoint-inverse (v1.0)
%   Author     : Ridwan Abdurrahman, Kakuta Laboratory, Department of Mechanical
%                Systems Engineering, Tokyo Metropolitan University
%   License    : MIT (see LICENSE)

clear; clc; close all;

%% === Settings ===============================================================
condition    = 'T28';         % data set, as in step 1 (T20, T28, T37)
start_dir    = pwd;
repo_root    = fileparts(mfilename('fullpath'));
if isempty(repo_root), repo_root = start_dir; end
res_dir      = fullfile(repo_root, 'results', condition);
assert(exist(fullfile(res_dir, 'step1_initial_estimates.mat'), 'file') == 2, ...
       'Run step1_initial_estimates.m with condition = %s first (no results in %s).', condition, res_dir);
cd(res_dir);
data_file    = 'step1_initial_estimates.mat';   % output of step 1
obs_source   = 'measured';    % 'measured' = NIR maps (C*_obs), 'fitted' = erfc/Gaussian fits
bc_source    = 'fitted';      % boundary values taken from the fitted fields
refit_profiles = true;        % true: least-squares fit of the profiles here (same method as step 1)
D0_source    = 'step1';       % 'step1' = D from Eqs. (20)-(21) in step 1, 'direct' = Eq. (20) here
sigma_S      = 0.03;          % mm, Gaussian width of the initial S (sigma* = 3)
cg_methods   = {'FR', 'PR', 'PB', 'DY'};
ref_method   = 'FR';          % reference CG method for part D, the sequential recovery and the field
                              % figures. DY reaches the lowest J, but it gave the least accurate fields
                              % in the synthetic test (Table 3), so FR is used. Set 'DY' to switch.
setappdata(0, 'ref_method', ref_method);   % read by the plotting functions (ref_method_index)
alpha_mode   = 'fixed';       % 'lcurve'  = one alpha for D and S from the joint L-curve
                              % 'lcurve2' = separate alphaD and alphaS: joint L-curve first, then
                              %             alternating L-curves (alphaS fixed -> alphaD, alphaD fixed -> alphaS)
                              % 'fixed'   = use alphaD/alphaS below (default: the values chosen
                              %             from the separate L-curves)
alpha_rule   = 'corner';      % 'corner' (maximum curvature) or 'discrepancy' (misfit = noise level)
alphas       = 10.^(-8:1:-3); % values tried for the joint L-curve
alphasD      = 10.^(-7:1:-2); % alphaD values for the separate L-curves (alpha_mode = 'lcurve2')
alphasS      = 10.^(-6:1:-1); % alphaS values for the separate L-curves
n_cycles_alpha = 2;          % alternating cycles (stops early when alphaD and alphaS no longer change)
noise_sd     = 0.087;         % mol/L, NIR concentration resolution (discrepancy principle)
alphaD       = 1e-5;          % used only when alpha_mode = 'fixed': lowest data misfit of the alphaD
                              % L-curve (Section 6.2); the same alphaD and alphaS are used for all
                              % data sets, because their L-curve corners differ (Section 6.6)
alphaS       = 1e-3;          % end of the steep branch of the alphaS L-curve; the curvature maximum
                              % (alphaS = 1e-2) makes S nearly uniform across the channel
line_search  = 'analytical';  % 'analytical' (+ backtracking, as in the article) or 'armijo'
Dref         = 1e-3;          % mm^2/s   (scaling only)
Sref         = 1e-2;          % mol/(L s) (scaling only)
Dmin         = 1e-5;          % mm^2/s, lower bound for D
Dmax         = 1e-1;          % mm^2/s, upper bound for D
run_gradient_check     = true;
run_sequential_control = true;
maxIter      = 300;           % CG iterations: 300 for all data sets, as in the article
                              % (T28: Sections 6.1-6.5; T20, T28, T37: Section 6.6)
cgopt = struct('maxIter', maxIter, 'tol_rel', 1e-4, 'ftol', 1e-7, 'stag_iter', 10, 'verbose', true, ...
               'line_search', line_search);

% --- sensitivity and uncertainty analyses (part D) ---
run_analyses   = true;        % false = stop after part C
lcurve_only    = false;       % true: run only the L-curves (alpha_mode 'lcurve' or 'lcurve2'), write
                              % step2_lcurve_results.txt and stop (results_summary.txt is not changed)
ana_method     = ref_method;  % CG method used for part D and the sequential recovery
sigmas         = [0.015 0.03 0.06];   % mm, sensitivity of the initial S shape
N_MC           = 20;          % Monte Carlo realisations
band_halfwidth = 0.1;         % mm, |x| < band_halfwidth is the "interfacial region"
run_items = struct('obs', true, 'sigma', true, 'mc', true);   % switch single analyses on/off

%% === A. Load step 1 results, fit profiles ====================================
L = load(data_file);
v_mean = L.v_dimensional;      % mm/s, mean velocity of this data set (from step 1)
D_lit  = L.D_lit_array;        % mm^2/s, reference values at the temperature of this data set
setappdata(0, 'D_lit', D_lit);
fprintf('Condition %s: v = %.2f mm/s, T = %d °C\n', L.condition, v_mean, L.T_C);
X = L.X;  Y = L.Y;                       % mm, Ny x Nx (meshgrid)
x = X(1, :);  y = Y(:, 1);
[Ny, Nx] = size(X);
dA = (x(2) - x(1))*(y(2) - y(1));
fprintf('Grid: Ny = %d (y), Nx = %d (x), dx = %.4f mm, dy = %.4f mm\n', Ny, Nx, x(2)-x(1), y(2)-y(1));
species = {'H2SO4', 'NaOH', 'Na2SO4'};
band = abs(X) < band_halfwidth;

if lcurve_only
    assert(any(strcmpi(alpha_mode, {'lcurve', 'lcurve2'})), 'lcurve_only needs alpha_mode = ''lcurve'' or ''lcurve2''.');
    fid = fopen('step2_lcurve_results.txt', 'w');
else
    fid = fopen('results_summary.txt', 'w');
end
out = @(varargin) fprintf_both(fid, varargin{:});
out('RESULTS OF STEP 2, condition %s  (%s)\n', condition, datestr(now));
out('line search %s, maxIter %d, alpha_mode %s, analysis method %s, D0 from %s\n\n', ...
    line_search, cgopt.maxIter, alpha_mode, ana_method, D0_source);

C_meas = {L.C1_obs, L.C2_obs, L.C3_obs};  % NIR maps on the computational grid (mol/L)
C_old  = {L.C1, L.C2, L.C3};              % profiles stored by step 1
C_fit  = C_old;
P_fit  = cell(1, 3);
if refit_profiles
    models = {'erfc_left', 'erfc_right', 'gauss'};
    fprintf('\nFitting erfc/Gaussian profiles to the measured maps (parameters linear in y)...\n');
    for s = 1:3
        [C_fit{s}, P_fit{s}] = fit_profiles(C_meas{s}, x, y, models{s});
    end
end
D_step1 = {L.D1_refined, L.D2_refined, L.D3_refined};   % mm^2/s, NaN outside the mask

if strcmpi(obs_source, 'measured'), C_obs = C_meas; else, C_obs = C_fit; end
if strcmpi(bc_source, 'fitted'),    C_bc  = C_fit;  else, C_bc  = C_meas; end
W = cellfun(@(c) double(~isnan(c)), C_obs, 'UniformOutput', false);

% --- quality of the profiles ---
out('--- quality of the profiles (vs measured maps, whole grid) ---\n');
out('%-7s | %-30s | %-30s | %s\n', 'species', 'profiles stored by step 1', 'least-squares fit (this run)', 'fit parameters (linear in y)');
for s = 1:3
    [r0, q0] = fit_stats(C_old{s}, C_meas{s});
    [r1, q1] = fit_stats(C_fit{s}, C_meas{s});
    if ~isempty(P_fit{s})
        ptxt = sprintf('a %.3f->%.3f, mu %.4f->%.4f mm, s %.4f->%.4f mm', P_fit{s}.a(1), P_fit{s}.a(end), ...
            P_fit{s}.mu(1), P_fit{s}.mu(end), P_fit{s}.s(1), P_fit{s}.s(end));
    else
        ptxt = '-';
    end
    out('%-7s | RMSE %.4f mol/L, R2 %6.3f   | RMSE %.4f mol/L, R2 %6.3f   | %s\n', species{s}, r0, q0, r1, q1, ptxt);
end
out('(NIR concentration resolution: %.3f mol/L)\n\n', noise_sd);
plot_fit_quality(x, y, C_meas, C_fit, C_old, species);

%% === A. Initial guess =========================================================
D0 = cell(1, 3);  S0 = cell(1, 3);
for s = 1:3
    S0{s} = init_source_balance(C_fit{s}, x, y, v_mean, D_step1{s}, sigma_S);
    if strcmpi(D0_source, 'direct')
        D0{s} = init_diffusion_direct(C_fit{s}, S0{s}, x, y, v_mean, Dmin, Dmax);
    else
        D0{s} = min(max(fill_nan_field(D_step1{s}, X, Y), Dmin), Dmax);
    end
    fprintf('Initial %-7s D: [%.3e, %.3e] mm^2/s,  S: [%.3e, %.3e] mol/(L s)\n', species{s}, ...
        min(D0{s}(:)), max(D0{s}(:)), min(S0{s}(:)), max(S0{s}(:)));
end

%% === A. Regularisation weights from the L-curves ==============================
opt = struct('alphaD', alphaD, 'alphaS', alphaS, 'Dref', Dref, 'Sref', Sref, 'Dmin', Dmin, 'Dmax', Dmax);
cgopt.Cpde = C_fit;          % forward error = residual of Eq. (2) with the smooth fitted fields
LC = [];  LC2 = [];
if any(strcmpi(alpha_mode, {'lcurve', 'lcurve2'}))
    fprintf('\nL-curve: %d inversions (alphaD = alphaS = alpha)...\n', numel(alphas));
    LC = lcurve_sweep(x, y, v_mean, C_obs, W, C_bc, opt, D0, S0, alphas, ana_method, cgopt, ...
                      noise_sd, alpha_rule, 'L_curve.png');
    alphaD = LC.alpha;  alphaS = LC.alpha;
    opt.alphaD = alphaD;  opt.alphaS = alphaS;
    out('--- L-curve (alphaD = alphaS = alpha) ---\n');
    for k = 1:size(LC.table, 1)
        out('  alpha %.1e: data misfit %.4e, H1 seminorm %.4e\n', LC.table(k, 1), LC.table(k, 2), LC.table(k, 3));
    end
    out('  L-curve corner (max curvature): alpha = %.1e\n', LC.alpha_corner);
    out('  discrepancy principle (misfit = %.4e at %.3f mol/L): alpha = %.1e\n', LC.target_misfit, noise_sd, LC.alpha_discrepancy);
    out('  selected (%s): alpha = %.1e\n\n', LC.rule, LC.alpha);
end
if strcmpi(alpha_mode, 'lcurve2')
    fprintf('\nSeparate alphaD and alphaS: alternating L-curves, starting from alpha = %.1e...\n', LC.alpha);
    LC2 = lcurve_sweep2(x, y, v_mean, C_obs, W, C_bc, opt, D0, S0, alphasD, alphasS, ana_method, cgopt, ...
                        LC.alpha, n_cycles_alpha, LC.runs, 'L_curve_separate.png');
    alphaD = LC2.alphaD;  alphaS = LC2.alphaS;
    opt.alphaD = alphaD;  opt.alphaS = alphaS;
    out('--- separate alphaD and alphaS (alternating L-curves) ---\n');
    for k = 1:numel(LC2.sweeps)
        sw = LC2.sweeps{k};
        if sw.which == 'D', fx = 'alphaS'; else, fx = 'alphaD'; end
        out('  cycle %d, alpha%s sweep (%s = %.1e fixed):\n', sw.cycle, sw.which, fx, sw.fixed);
        for j = 1:size(sw.table, 1)
            out('    alpha%s %.1e: data misfit %.4e, H1 seminorm of %s %.4e\n', sw.which, sw.table(j, 1), ...
                sw.table(j, 2), sw.which, sw.table(j, 3));
        end
        out('    -> corner: alpha%s = %.1e\n', sw.which, sw.selected);
    end
    out('  selected: alphaD = %.1e, alphaS = %.1e (%d cycle(s), %d inversions)\n\n', alphaD, alphaS, ...
        LC2.cycles, size(LC2.runs, 1));
end
out('Regularisation used for all runs: alphaD = %.1e, alphaS = %.1e (%s)\n\n', alphaD, alphaS, alpha_mode);
if lcurve_only
    save('step2_lcurve_results.mat', 'LC', 'LC2', 'alphas', 'alphasD', 'alphasS');
    fclose(fid);
    fprintf('\nlcurve_only = true: L-curves written to step2_lcurve_results.txt / step2_lcurve_results.mat; stopping here.\n');
    cd(start_dir);
    return
end
prob = make_problem(x, y, v_mean, C_obs, W, C_bc, opt);
mkprob = @(Cobs) make_problem(x, y, v_mean, Cobs, W, C_bc, opt);   % same alpha, other data

if run_gradient_check
    fprintf('\n');
    gc = gradient_check(D0, S0, prob, 2);
    out('Gradient check (FD vs adjoint), max relative error: D %.1e, S %.1e\n\n', max(gc(:, 1)), max(gc(:, 2)));
end

%% === B. Run all CG methods (simultaneous estimation) ==========================
all_results = cell(1, numel(cg_methods));
for m = 1:numel(cg_methods)
    fprintf('\n=============== %s ===============\n', cg_methods{m});
    all_results{m} = cg_optimize(prob, D0, S0, cg_methods{m}, cgopt);
end
out('--- CG methods (Table 1) ---\n');
out('%-6s %-20s %7s %12s %12s %8s %12s %12s %9s %8s\n', 'Method', 'Stop reason', 'Updates', ...
    'J initial', 'J final', 'Drop %', 'Grad norm', 'Fwd error', 'Restarts', 'Time s');
for m = 1:numel(cg_methods)
    r = all_results{m};
    out('%-6s %-20s %7d %12.5e %12.5e %8.2f %12.4e %12.4e %9d %8.1f\n', cg_methods{m}, r.stop_reason, ...
        r.num_updates, r.initial_cost, r.final_cost, 100*(1 - r.final_cost/r.initial_cost), ...
        r.final_grad_norm, r.forward_error_history(end), r.cg_restart_count, r.total_time);
end
out('\n');

% --- choice of the reference method (the lowest J is not the criterion) ---
out('--- Reference method for the fields, figures and analyses: %s (set by ref_method) ---\n', ref_method);
out('%-6s %12s %12s %12s %10s %12s %9s\n', 'Method', 'J final', 'data misfit', 'reg. term', 'rel. grad', 'fwd error', 'restarts');
for m = 1:numel(cg_methods)
    r = all_results{m};  Jd = r.J_data_history(end);
    out('%-6s %12.4e %12.4e %12.4e %10.3f %12.4e %9d\n', cg_methods{m}, r.final_cost, Jd, r.final_cost - Jd, ...
        r.final_grad_norm/r.grad_norm_history(1), r.forward_error_history(end), r.cg_restart_count);
end
[~, iJ] = min(cellfun(@(r) r.final_cost, all_results));
[~, iD] = min(cellfun(@(r) r.J_data_history(end), all_results));
out('lowest J: %s; lowest data misfit: %s. A lower J does not imply more accurate fields (synthetic test),\n', ...
    cg_methods{iJ}, cg_methods{iD});
out('so the reference method is set by ref_method and not chosen by the lowest J.\n\n');

% baseline for part D = the run of ana_method above (no extra run)
ib = find(strcmpi(cg_methods, ana_method), 1);
if isempty(ib)
    base = cg_optimize(prob, D0, S0, ana_method, cgopt);
else
    base = all_results{ib};
end
res_rms = sqrt(2*base.J_data_history(end)/(sum(cellfun(@(w) sum(w(:)), W))*dA));
out('--- Baseline (%s, measured maps, sigma = %.3f mm) ---\n', ana_method, sigma_S);
out('data misfit %.4e, RMS data residual %.4f mol/L; misfit expected at the NIR noise level %.4e\n', ...
    base.J_data_history(end), res_rms, 0.5*sum(cellfun(@(w) sum(w(:)), W))*noise_sd^2*dA);
% outside the reaction zone on the side where the species is present (as in Table 2)
bulk_mask = {X < -0.25, X > 0.25, abs(X) > 0.25};
bulk_lab  = {'x<-0.25', 'x>0.25', '|x|>0.25'};
for s = 1:3
    Dm = base.(sprintf('D%d', s));  Sm = base.(sprintf('S%d', s));
    out('  %-7s D bulk (%s mm) median %.3e, interfacial max %.3e mm2/s (lit. %.2e) | S min %.3e, max %.3e mol/(L s)\n', ...
        species{s}, bulk_lab{s}, median(Dm(bulk_mask{s})), max(Dm(band)), D_lit(s), min(Sm(:)), max(Sm(:)));
end
out('\n');

%% === B. Sequential recovery: D fixed at literature bulk values =========
if run_sequential_control
    Dlit = {D_lit(1)*ones(Ny, Nx), D_lit(2)*ones(Ny, Nx), D_lit(3)*ones(Ny, Nx)};
    cg_seq = cgopt;  cg_seq.fixD = true;
    fprintf('\n=============== Sequential recovery (D fixed, %s) ===============\n', ana_method);
    seq_result = cg_optimize(prob, Dlit, S0, ana_method, cg_seq);
    out('--- sequential recovery (D fixed at literature values, %s) ---\n', ana_method);
    out('  data misfit: simultaneous %.4e, sequential %.4e (ratio %.2f)\n', base.J_data_history(end), ...
        seq_result.J_data_history(end), seq_result.J_data_history(end)/base.J_data_history(end));
    for s = 1:3
        out('  %-7s peak |S|: simultaneous %.3e, sequential %.3e mol/(L s); dS %.1f%%\n', species{s}, ...
            max(abs(base.(sprintf('S%d', s))(:))), max(abs(seq_result.(sprintf('S%d', s))(:))), ...
            100*reldiff(seq_result, base, 'S', s, band));
    end
    out('\n');
end

%% === C. Figures and tables ========
create_method_comparison_plots(all_results, cg_methods);
create_line_profile_plots(all_results, cg_methods, X, Y);
create_parameter_comparison_plots_physical(all_results, cg_methods, X, Y);
create_performance_comparison_table(all_results, cg_methods);
plot_change_maps(all_results, cg_methods, D0, S0, X, Y);
if run_sequential_control
    plot_sequential_comparison(all_results, cg_methods, seq_result, X, Y);
end

save('step2_cg_results.mat', 'all_results', 'cg_methods', 'D0', 'S0', 'X', 'Y', 'x', 'y', 'LC', 'LC2', ...
     'alpha_mode', 'line_search', 'alphaD', 'alphaS', 'Dref', 'Sref', 'obs_source', 'v_mean', 'sigma_S', 'C_fit');
if run_sequential_control, save('step2_cg_results.mat', 'seq_result', '-append'); end
fprintf('\nSaved step2_cg_results.mat\n');

%% === D. Sensitivity and uncertainty analyses ===============================
if run_analyses
    % --- fitted fields instead of measured maps in the misfit ---
    if run_items.obs
        fprintf('\nRun with fitted fields in the misfit...\n');
        rf = cg_optimize(mkprob(C_fit), D0, S0, ana_method, quiet(cgopt));
        out('--- misfit on fitted fields vs measured maps (relative difference in |x| < %.2f mm) ---\n', band_halfwidth);
        for s = 1:3
            out('  %-7s dD %.1f%%   dS %.1f%%\n', species{s}, 100*reldiff(rf, base, 'D', s, band), 100*reldiff(rf, base, 'S', s, band));
        end
        out('\n');
    end

    % --- Gaussian width of the initial S ---
    if run_items.sigma
        out('--- sensitivity to sigma of the initial S (relative to sigma = %.3f mm) ---\n', sigma_S);
        for sg = sigmas
            if abs(sg - sigma_S) < 1e-12, continue; end
            fprintf('Run with sigma = %.3f mm...\n', sg);
            S0s = cell(1, 3);
            for s = 1:3, S0s{s} = init_source_balance(C_fit{s}, x, y, v_mean, D_step1{s}, sg); end
            rs = cg_optimize(prob, D0, S0s, ana_method, quiet(cgopt));
            line = sprintf('  sigma %.3f mm: J_end %.4e |', sg, rs.final_cost);
            for s = 1:3
                Ia = trapz(x, rs.(sprintf('S%d', s)), 2);  Ibs = trapz(x, base.(sprintf('S%d', s)), 2);
                line = [line, sprintf(' %s dS %.1f%%, d(int S dx) %.1f%% |', species{s}, ...
                    100*reldiff(rs, base, 'S', s, band), 100*norm(Ia - Ibs)/norm(Ibs))]; %#ok<AGROW>
            end
            out('%s\n', line);
        end
        out('\n');
    end

    % --- Monte Carlo uncertainty (parametric bootstrap) ---
    if run_items.mc
        out('--- Monte Carlo, %d realisations, noise %.4f mol/L added to the baseline model fields ---\n', N_MC, res_rms);
        Ds = zeros(Ny, Nx, 3, N_MC);  Ss = Ds;
        for n = 1:N_MC
            fprintf('Monte Carlo %d/%d...\n', n, N_MC);
            Cb = cellfun(@(c) c + res_rms*randn(size(c)), base.C_model, 'UniformOutput', false);
            r = cg_optimize(mkprob(Cb), D0, S0, ana_method, quiet(cgopt));
            for s = 1:3
                Ds(:, :, s, n) = r.(sprintf('D%d', s));  Ss(:, :, s, n) = r.(sprintf('S%d', s));
            end
        end
        figH = figure('Color', 'w', 'Position', [60 60 1300 600]);
        for s = 1:3
            sdD = std(Ds(:, :, s, :), 0, 4);  mD = mean(Ds(:, :, s, :), 4);
            sdS = std(Ss(:, :, s, :), 0, 4);  mS = mean(Ss(:, :, s, :), 4);
            relD = median(sdD(band)./abs(mD(band)));
            sSc = max(abs(mS(band)));  big = band & abs(mS) > 0.2*sSc;
            relS = median(sdS(big)./abs(mS(big)));
            out('  %-7s median std/|mean| in |x| < %.2f mm: D %.1f%%, S %.1f%% (S where |S| > 20%% of its peak)\n', ...
                species{s}, band_halfwidth, 100*relD, 100*relS);
            subplot(2, 3, s); imagesc(x, y, sdD); axis xy; colorbar; title(sprintf('std D %s (mm^2/s)', species{s}));
            subplot(2, 3, 3 + s); imagesc(x, y, sdS); axis xy; colorbar; title(sprintf('std S %s (mol/(L s))', species{s}));
        end
        save_figure_safe(figH, 'MC_std.png', 300);
        save('step2_cg_results.mat', 'Ds', 'Ss', '-append');
        out('\n');
    end
end
fclose(fid);
fprintf('\nAll numbers written to results_summary.txt\n');

%% === Optional: display final recommendation ==============
display_method_recommendation(all_results, cg_methods);
fprintf('Results of %s are in %s\n', condition, res_dir);
cd(start_dir);

%% ============================================================================
%% LOCAL FUNCTIONS: SOLVER, ADJOINT, PROFILE FIT, L-CURVE AND CG
%% ============================================================================

function prob = make_problem(x, y, v, Cobs, W, Cbc, opt)
% MAKE_PROBLEM  Collect grid, data and settings for the inverse problem.
%   x (1 x Nx, mm), y (Ny x 1 or 1 x Ny, mm), v (mm/s, scalar or 1 x Nx)
%   Cobs{s}  observed concentration used in the data misfit (mol/L), Ny x Nx
%   W{s}     misfit weights (1 = use point, 0 = ignore, e.g. NaN in data)
%   Cbc{s}   field whose inlet row and wall columns give the Dirichlet values
%   opt      struct: alphaD, alphaS, Dref, Sref, Dmin, Dmax

    prob.x = x(:).';  prob.y = y(:);
    prob.Nx = numel(x);  prob.Ny = numel(y);
    prob.dx = x(2) - x(1);  prob.dy = y(2) - y(1);
    prob.v  = v;
    prob.nsp = numel(Cobs);

    bnd = false(prob.Ny, prob.Nx);
    bnd(1, :) = true;          % inlet
    bnd(:, 1) = true;          % wall x = x_min
    bnd(:, end) = true;        % wall x = x_max
    prob.bnd = bnd;
    prob.interior = ~bnd;

    for s = 1:prob.nsp
        c = Cobs{s};
        w = double(W{s});
        w(isnan(c)) = 0;
        c(isnan(c)) = 0;
        prob.Cobs{s} = c;
        prob.W{s} = w;
        g = Cbc{s};
        g(isnan(g)) = 0;
        prob.gbc{s} = g(bnd);   % Dirichlet values on boundary nodes
    end

    prob.alphaD = opt.alphaD;  prob.alphaS = opt.alphaS;
    prob.Dref = opt.Dref;      prob.Sref = opt.Sref;
    prob.Dmin = opt.Dmin;      prob.Dmax = opt.Dmax;
    prob.K = reg_matrix(prob.Ny, prob.Nx, prob.dx, prob.dy);
end

function A = adr_assemble(D, prob)
% ADR_ASSEMBLE  Sparse matrix of the discrete Eq. (2) for one species:
%
%     v(x) dC/dy = d/dx( D dC/dx ) + S          (steady, depth-averaged)
%
% Discretisation (conservative, finite volume in x, implicit march in y):
%   v_j (C_ij - C_{i-1,j})/dy
%     - [ D_{i,j+1/2}(C_{i,j+1}-C_ij) - D_{i,j-1/2}(C_ij-C_{i,j-1}) ]/dx^2 = S_ij
%   with D_{j+1/2} = (D_j + D_{j+1})/2.
% Rows for boundary nodes (inlet y = 0 and the two walls x = x_min, x_max)
% are identity rows: C = C_bc (Dirichlet values taken from the data).
% Outlet (y = y_max) needs no condition because the problem is parabolic in y.
%
% Indexing: MATLAB column-major, k = i + (j-1)*Ny, i = y-index, j = x-index.

    Ny = prob.Ny;  Nx = prob.Nx;  dx = prob.dx;  dy = prob.dy;
    N  = Ny*Nx;
    idx = reshape(1:N, Ny, Nx);

    v = prob.v;                        % scalar or 1 x Nx
    if isscalar(v), v = v*ones(1, Nx); end
    V = repmat(v(:).', Ny, 1);

    Dp = 0.5*(D(:, 2:Nx-1) + D(:, 3:Nx));     % D_{j+1/2}, interior columns
    Dm = 0.5*(D(:, 2:Nx-1) + D(:, 1:Nx-2));   % D_{j-1/2}
    Dp = Dp(2:Ny, :);  Dm = Dm(2:Ny, :);
    Vi = V(2:Ny, 2:Nx-1);

    kc = idx(2:Ny,   2:Nx-1);   % centre
    ks = idx(1:Ny-1, 2:Nx-1);   % upstream (i-1)
    ke = idx(2:Ny,   3:Nx);     % j+1
    kw = idx(2:Ny,   1:Nx-2);   % j-1

    rows = [kc(:); kc(:); kc(:); kc(:)];
    cols = [kc(:); ks(:); ke(:); kw(:)];
    vals = [Vi(:)/dy + (Dp(:) + Dm(:))/dx^2;
            -Vi(:)/dy;
            -Dp(:)/dx^2;
            -Dm(:)/dx^2];

    kb = idx(prob.bnd);                       % boundary nodes -> identity rows
    rows = [rows; kb(:)];  cols = [cols; kb(:)];  vals = [vals; ones(numel(kb), 1)];

    A = sparse(rows, cols, vals, N, N);
end

function [C, A] = forward_solve(D, S, s, prob)
% FORWARD_SOLVE  Solve the discrete Eq. (2) for species s.
    A = adr_assemble(D, prob);
    b = S(:);
    b(prob.bnd(:)) = prob.gbc{s};
    C = reshape(A \ b, prob.Ny, prob.Nx);
end

function [J, gD, gS, info] = objective_gradient(D, S, prob, need_grad)
% OBJECTIVE_GRADIENT  Eq. (3) and its exact (discrete-adjoint) gradient.
%
%   J = 1/2 sum_s sum W (C_s - Cobs_s)^2 dx dy
%     + alphaD/2 sum_s |grad(D_s/Dref)|^2 dx dy + alphaS/2 sum_s |grad(S_s/Sref)|^2 dx dy
%
%   D, S : cell arrays {1..nsp} of Ny x Nx fields (D in mm^2/s, S in mol/(L s))
%   gD, gS : L2 gradients (per unit area) w.r.t. D and S, same units as dJ/dD, dJ/dS
%
% Adjoint: A(D)' * lambda = W (C - Cobs) dx dy   (exact transpose of the forward matrix,
% so it contains the reversed advection term automatically).
    if nargin < 4, need_grad = true; end
    dx = prob.dx;  dy = prob.dy;  dA = dx*dy;
    Ny = prob.Ny;  Nx = prob.Nx;
    Jd = 0;  Jr = 0;
    gD = cell(1, prob.nsp);  gS = cell(1, prob.nsp);  C = cell(1, prob.nsp);
    for s = 1:prob.nsp
        [Cs, A] = forward_solve(D{s}, S{s}, s, prob);
        C{s} = Cs;
        r = prob.W{s} .* (Cs - prob.Cobs{s});
        Jd = Jd + 0.5*sum(r(:) .* (Cs(:) - prob.Cobs{s}(:)))*dA;

        tD = D{s}(:)/prob.Dref;  tS = S{s}(:)/prob.Sref;
        KD = prob.K*tD;  KS = prob.K*tS;
        Jr = Jr + 0.5*prob.alphaD*(tD.'*KD)*dA + 0.5*prob.alphaS*(tS.'*KS)*dA;

        if need_grad
            lam = reshape(A.' \ (r(:)*dA), Ny, Nx);
            lam(prob.bnd) = 0;                    % boundary rows do not depend on D or S
            % dJ/dS (nodal) = lambda on interior nodes
            gS_nodal = lam;
            % dJ/dD (nodal) from the flux terms, faces f = (j, j+1):
            % dJ/dD_m = [ (lam_m - lam_{m+1}) E_m + (lam_{m-1} - lam_m) E_{m-1} ] / (2 dx^2)
            E  = Cs(:, 2:Nx) - Cs(:, 1:Nx-1);            % Ny x (Nx-1)
            dl = lam(:, 1:Nx-1) - lam(:, 2:Nx);          % lam_j - lam_{j+1}
            P  = dl .* E / (2*dx^2);                     % face contribution
            gD_nodal = zeros(Ny, Nx);
            gD_nodal(:, 1:Nx-1) = gD_nodal(:, 1:Nx-1) + P;   % face (m, m+1) -> node m
            gD_nodal(:, 2:Nx)   = gD_nodal(:, 2:Nx)   + P;   % face (m-1, m) -> node m
            gD_nodal(1, :) = 0;   % inlet row: D there is not used by any interior row
            % convert to L2 gradients and add regularisation
            gD{s} = gD_nodal/dA + reshape(prob.alphaD*KD, Ny, Nx)/prob.Dref;
            gS{s} = gS_nodal/dA + reshape(prob.alphaS*KS, Ny, Nx)/prob.Sref;
        end
    end
    J = Jd + Jr;
    info.Jdata = Jd;  info.Jreg = Jr;  info.C = C;
end

function K = reg_matrix(Ny, Nx, dx, dy)
% REG_MATRIX  K = Gx'Gx + Gy'Gy with forward-difference gradients (column-major).
% Then  0.5*sum(|grad f|^2)*dx*dy = 0.5*f(:)'*K*f(:)*dx*dy   (H1 seminorm)
% and its L2 gradient is K*f(:)  (a Neumann Laplacian with a minus sign).
    ex = ones(Nx, 1);  ey = ones(Ny, 1);
    Dx1 = spdiags([-ex ex], [0 1], Nx-1, Nx) / dx;   % (Nx-1) x Nx
    Dy1 = spdiags([-ey ey], [0 1], Ny-1, Ny) / dy;   % (Ny-1) x Ny
    Gx = kron(Dx1, speye(Ny));     % derivative along x (slow index)
    Gy = kron(speye(Nx), Dy1);     % derivative along y (fast index)
    K = Gx.'*Gx + Gy.'*Gy;
end

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

function out = gradient_check(D, S, prob, ntests)
% GRADIENT_CHECK  Compare the adjoint gradient with central finite differences
% along random smooth directions. Relative errors ~1e-6 or smaller mean the
% gradient is exact for the discrete problem.
    if nargin < 4, ntests = 3; end
    dA = prob.dx*prob.dy;
    [~, gD, gS] = objective_gradient(D, S, prob, true);
    out = zeros(ntests, 2);
    for t = 1:ntests
        for which = 1:2            % 1 = perturb D, 2 = perturb S
            pD = cell(1, prob.nsp);  pS = cell(1, prob.nsp);
            dirdot = 0;
            for s = 1:prob.nsp
                R = smooth_random(prob.Ny, prob.Nx);
                if which == 1
                    pD{s} = R*prob.Dref;  pS{s} = zeros(prob.Ny, prob.Nx);
                    dirdot = dirdot + sum(gD{s}(:).*pD{s}(:))*dA;
                else
                    pD{s} = zeros(prob.Ny, prob.Nx);  pS{s} = R*prob.Sref;
                    dirdot = dirdot + sum(gS{s}(:).*pS{s}(:))*dA;
                end
            end
            h = 1e-4;
            Dp = cellfun(@(a, b) a + h*b, D, pD, 'UniformOutput', false);
            Dm = cellfun(@(a, b) a - h*b, D, pD, 'UniformOutput', false);
            Sp = cellfun(@(a, b) a + h*b, S, pS, 'UniformOutput', false);
            Sm = cellfun(@(a, b) a - h*b, S, pS, 'UniformOutput', false);
            fd = (objective_gradient(Dp, Sp, prob, false) - objective_gradient(Dm, Sm, prob, false))/(2*h);
            out(t, which) = abs(fd - dirdot)/max(abs(fd), eps);
        end
    end
    fprintf('Gradient check (relative error, FD vs adjoint):\n');
    fprintf('   D directions: %s\n', sprintf('%.2e  ', out(:, 1)));
    fprintf('   S directions: %s\n', sprintf('%.2e  ', out(:, 2)));
end

function R = smooth_random(Ny, Nx)
    R = randn(Ny, Nx);
    k = ones(5)/25;
    R = conv2(R, k, 'same');
    R = R/max(abs(R(:)));
end

function S0 = init_source_balance(C, x, y, v, D0, sigma, x0)
% INIT_SOURCE_BALANCE  Initial source term from the width-integrated Eq. (2),
% Eqs. (17)-(19), in physical units.
%   Phi(y) = d/dy int v C dx - [ D dC/dx ]_{x1}^{x2}        (mol mm /(L s))
%   S0(x,y) = Phi(y) * G(x),  G = Gaussian with width sigma (mm), int G dx = 1
%   C: Ny x Nx fitted concentration (mol/L), x (mm), y (mm), v (mm/s),
%   D0: Ny x Nx diffusion field used for the wall flux (mm^2/s, NaN allowed),
%   sigma: Gaussian width in mm (sigma* = 3 in the article = 0.03 mm), x0: centre.
    if nargin < 7, x0 = 0; end
    x = x(:).';  y = y(:);
    Nx = numel(x);  dx = x(2) - x(1);
    if isscalar(v), v = v*ones(1, Nx); end
    Iv = trapz(x, C.*repmat(v(:).', size(C, 1), 1), 2);
    Iv = smooth_ma(Iv, 5);
    dIv = gradient(Iv, y(2) - y(1));
    dCdx_L = (-3*C(:, 1) + 4*C(:, 2) - C(:, 3))/(2*dx);
    dCdx_R = ( 3*C(:, end) - 4*C(:, end-1) + C(:, end-2))/(2*dx);
    Dfill = median(D0(~isnan(D0)));
    if isempty(Dfill) || isnan(Dfill), Dfill = 1e-3; end     % mm^2/s fallback
    DL = D0(:, 1);  DR = D0(:, end);
    DL(isnan(DL)) = Dfill;  DR(isnan(DR)) = Dfill;
    Phi = dIv - (DR.*dCdx_R - DL.*dCdx_L);
    G = exp(-(x - x0).^2/(2*sigma^2));
    G = G/trapz(x, G);
    S0 = Phi*G;
end

function z = smooth_ma(z, n)
    k = ones(n, 1)/n;
    zp = [repmat(z(1), floor(n/2), 1); z; repmat(z(end), floor(n/2), 1)];
    z = conv(zp, k, 'valid');
end

function D0 = init_diffusion_direct(C, S, x, y, v, Dmin, Dmax, thr)
% INIT_DIFFUSION_DIRECT  Direct pointwise inversion of Eq. (2) (Eq. (20)):
%     D = ( v dC/dy - S ) / ( d2C/dx2 )
% evaluated with the smooth fitted field C (mol/L), x and y in mm, v in mm/s,
% S in mol/(L s). Points where |d2C/dx2| < thr*max|d2C/dx2| or where D falls
% outside [Dmin, Dmax] are set to NaN and then filled by interpolation.
    if nargin < 8, thr = 0.05; end
    x = x(:).';  y = y(:);
    dx = x(2) - x(1);  dy = y(2) - y(1);
    [dCdx, dCdy] = gradient(C, dx, dy); %#ok<ASGLU>
    d2 = zeros(size(C));
    d2(:, 2:end-1) = (C(:, 3:end) - 2*C(:, 2:end-1) + C(:, 1:end-2))/dx^2;
    if isscalar(v), V = v; else, V = repmat(v(:).', size(C, 1), 1); end
    D0 = (V.*dCdy - S)./d2;
    bad = abs(d2) < thr*max(abs(d2(:))) | ~isfinite(D0) | D0 < Dmin | D0 > Dmax;
    bad(:, [1 end]) = true;
    D0(bad) = NaN;
    if all(isnan(D0(:)))
        D0(:) = 1e-3;
    else
        [X, Y] = meshgrid(x, y);
        D0 = fill_nan_field(D0, X, Y);
    end
    D0 = min(max(D0, Dmin), Dmax);
end

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

function [a_corner, a_disc, at_edge] = lcurve_choice(Lc, target)
% corner = maximum curvature of (log misfit, log seminorm) vs log alpha, taken over the
% INTERIOR points only: at the two end points the curvature comes from one-sided differences
% and is not reliable, so an end point is never returned as the corner.
% at_edge = true when an end point has a larger (unreliable) curvature than the interior
% corner, i.e. the range should be extended to confirm the corner.
% discrepancy = alpha whose data misfit is closest to the noise-level target
    la = log10(Lc(:, 1));  lr = log10(Lc(:, 2));  ls = log10(max(Lc(:, 3), realmin));
    at_edge = false;
    if numel(la) >= 3
        r1 = gradient(lr, la);  s1 = gradient(ls, la);
        r2 = gradient(r1, la);  s2 = gradient(s1, la);
        kap = (r1.*s2 - r2.*s1)./max((r1.^2 + s1.^2).^1.5, eps);
        [kmax, i] = max(kap(2:end-1));
        a_corner = Lc(i + 1, 1);
        at_edge = max(kap([1, end])) > kmax;
    else
        a_corner = NaN;
    end
    [~, j] = min(abs(log10(Lc(:, 2)) - log10(target)));
    a_disc = Lc(j, 1);
end

function lc_save_png(figH, fname)
    if exist('exportgraphics', 'file') == 2
        exportgraphics(figH, fname, 'Resolution', 300);
    else
        print(figH, fname, '-dpng', '-r200');
    end
end

function [etaD, etaS] = scaled_seminorms(r, prob)
% SCALED_SEMINORMS  H1 seminorms of the scaled fields of a CG result r:
%   etaD = 0.5*sum_i |grad(D_i/Dref)|^2 dA,  etaS = 0.5*sum_i |grad(S_i/Sref)|^2 dA,
%   so that J_reg = alphaD*etaD + alphaS*etaS (same definition as objective_gradient).
    dA = prob.dx*prob.dy;  etaD = 0;  etaS = 0;
    for s = 1:prob.nsp
        tD = r.(sprintf('D%d', s))(:)/prob.Dref;
        tS = r.(sprintf('S%d', s))(:)/prob.Sref;
        etaD = etaD + 0.5*(tD.'*(prob.K*tD))*dA;
        etaS = etaS + 0.5*(tS.'*(prob.K*tS))*dA;
    end
end

function F = fill_nan_field(F, X, Y)
% FILL_NAN_FIELD  Fill NaN points by linear interpolation from valid points,
% nearest-neighbour outside their convex hull (same idea as Section 5).
    bad = isnan(F);
    if ~any(bad(:)), return; end
    good = ~bad;
    F(bad) = griddata(X(good), Y(good), F(good), X(bad), Y(bad), 'linear');
    bad2 = isnan(F);
    if any(bad2(:))
        good2 = ~bad2;
        F(bad2) = griddata(X(good2), Y(good2), F(good2), X(bad2), Y(bad2), 'nearest');
    end
end

%% ============================================================================
%% LOCAL FUNCTIONS: FIGURES AND HELPERS FOR THE SENSITIVITY ANALYSES
%% ============================================================================

function plot_change_maps(all_results, cg_methods, D0, S0, X, Y)
% Change of D and S produced by the optimisation (final minus initial),
% for the reference method (ref_method in the main script).
    b = ref_method_index(all_results, cg_methods);
    R = all_results{b};
    names = {'H_2SO_4', 'NaOH', 'Na_2SO_4'};
    figH = figure('Name', sprintf('Change by optimisation (%s)', cg_methods{b}), ...
                  'Position', [60, 60, 1400, 650], 'Color', 'w');
    for s = 1:3
        subplot(2, 3, s);
        imagesc(X(1, :), Y(:, 1), R.(sprintf('D%d', s)) - D0{s});
        axis xy; colorbar; xlabel('x (mm)'); ylabel('y (mm)');
        title(sprintf('\\DeltaD %s (mm^2/s)', names{s}));
        subplot(2, 3, 3 + s);
        imagesc(X(1, :), Y(:, 1), R.(sprintf('S%d', s)) - S0{s});
        axis xy; colorbar; xlabel('x (mm)'); ylabel('y (mm)');
        title(sprintf('\\DeltaS %s (mol/(L s))', names{s}));
    end
    save_figure_safe(figH, 'change_by_optimisation.png', 300);
end

function plot_sequential_comparison(all_results, cg_methods, seq, X, Y)
% Simultaneous (reference CG method) vs sequential recovery (D fixed at bulk value):
% S line profiles at y = 1, 4, 8 mm.
    b = ref_method_index(all_results, cg_methods);
    R = all_results{b};
    x = X(1, :);  yv = Y(:, 1);
    names = {'H_2SO_4', 'NaOH', 'Na_2SO_4'};
    yt = [1 4 8];  ls = {'-', '--', ':'};
    figH = figure('Name', 'Simultaneous vs sequential', 'Position', [60, 60, 1400, 380], 'Color', 'w');
    for s = 1:3
        subplot(1, 3, s); hold on;
        for k = 1:3
            [~, i] = min(abs(yv - yt(k)));
            plot(x, R.(sprintf('S%d', s))(i, :), ['b' ls{k}], 'LineWidth', 1.5);
            plot(x, seq.(sprintf('S%d', s))(i, :), ['r' ls{k}], 'LineWidth', 1.5);
        end
        hold off; grid on; xlabel('x (mm)'); ylabel('S (mol/(L s))'); title(names{s});
        if s == 1
            legend({'simultaneous y=1', 'sequential y=1', 'simultaneous y=4', 'sequential y=4', ...
                    'simultaneous y=8', 'sequential y=8'}, 'Location', 'best');
        end
    end
    fprintf('Data misfit: simultaneous (%s) %.5e, sequential %.5e\n', cg_methods{b}, ...
        R.J_data_history(end), seq.J_data_history(end));
    save_figure_safe(figH, 'simultaneous_vs_sequential_S.png', 300);
end

function [idx, name] = ref_method_index(all_results, cg_methods)
% Index of the reference CG method used for the field figures. It is the method set by
% ref_method in the main script (stored with setappdata). Only if none is set, the
% method with the lowest final J is used.
    name = '';
    if isappdata(0, 'ref_method'), name = getappdata(0, 'ref_method'); end
    idx = find(strcmpi(cg_methods, name), 1);
    if isempty(idx)
        [~, idx] = min(cellfun(@(r) r.final_cost, all_results));
        warning('ref_method not set or not run; using the method with the lowest J (%s).', cg_methods{idx});
    end
    name = cg_methods{idx};
end

function fprintf_both(fid, varargin)
% print to the command window and to the results file
    fprintf(varargin{:});
    fprintf(fid, varargin{:});
end

function c = quiet(c)
% same CG options without the per-iteration printout
    c.verbose = false;
end

function [rmse, R2] = fit_stats(Cf, Cm)
    m = ~isnan(Cm) & ~isnan(Cf);
    r = Cf(m) - Cm(m);
    rmse = sqrt(mean(r.^2));
    R2 = 1 - sum(r.^2)/sum((Cm(m) - mean(Cm(m))).^2);
end

function d = reldiff(ra, rb, field, s, band)
% relative L2 difference of field (D or S) of species s inside the band
    A = ra.(sprintf('%s%d', field, s));  B = rb.(sprintf('%s%d', field, s));
    d = norm(A(band) - B(band))/max(norm(B(band)), eps);
end

function plot_fit_quality(x, y, C_meas, C_fit, C_old, species)
% measured vs fitted transverse profiles at y = 1, 4, 8 mm (also checks the x alignment)
    figH = figure('Color', 'w', 'Position', [60 60 1400 380]);
    yt = [1 4 8];  ls = {'-', '--', ':'};  cols = {[0.85 0.1 0.1], [0.1 0.3 0.85], [0.1 0.6 0.2]};
    for s = 1:3
        subplot(1, 3, s); hold on;
        for k = 1:3
            [~, i] = min(abs(y - yt(k)));
            plot(x, C_meas{s}(i, :), '.', 'Color', cols{s}*0.6 + 0.4, 'MarkerSize', 6);
            plot(x, C_fit{s}(i, :), ls{k}, 'Color', cols{s}, 'LineWidth', 1.6);
            plot(x, C_old{s}(i, :), ls{k}, 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8);
        end
        hold off; grid on; xlabel('x (mm)'); ylabel('C (mol/L)'); title(species{s});
    end
    annotation('textbox', [0.01 0.92 0.98 0.08], 'String', ...
        'dots: measured;  coloured lines: least-squares fit;  grey lines: profiles stored by step 1;  solid/dashed/dotted: y = 1/4/8 mm', ...
        'EdgeColor', 'none', 'FontSize', 8);
    save_figure_safe(figH, 'fit_quality.png', 300);
end

%% ============================================================================
%% LOCAL FUNCTIONS: FIGURE EXPORT AND PLOTTING
%% ============================================================================

function save_figure_safe(figHandle, fname, dpi)
    % Figure export: exportgraphics, then print with painters, then a low-resolution
    % print as a last resort. Never throws an error.
    if nargin < 3, dpi = 600; end
    if isempty(figHandle) || ~isgraphics(figHandle)
        warning('save_figure_safe: invalid figure handle for %s', fname);
        return;
    end

    drawnow;  % finish rendering before export

    ok = false;

    % --- exportgraphics (MATLAB R2020a or newer)
    if exist('exportgraphics', 'file') == 2
        try
            exportgraphics(figHandle, fname, 'Resolution', dpi);
            ok = true;
        catch ME
            fprintf('   exportgraphics failed (%s); trying print...\n', ME.message);
        end
    end

    % --- Fallback 1: print with the painters renderer
    if ~ok
        try
            set(figHandle, 'Renderer', 'painters');
            print(figHandle, fname, '-dpng', sprintf('-r%d', dpi));
            ok = true;
        catch ME
            fprintf('   print(painters) failed (%s); trying low-DPI...\n', ME.message);
        end
    end

    % --- Fallback 2: OpenGL at low resolution
    if ~ok
        try
            set(figHandle, 'Renderer', 'opengl');
            print(figHandle, fname, '-dpng', '-r150');
            ok = true;
        catch ME2
            fprintf('   ALL export methods failed for %s: %s\n', fname, ME2.message);
        end
    end

    if ok
        fprintf('   Exported -> %s\n', fname);
    end
end

function create_method_comparison_plots(all_results, cg_methods)
    num_methods = length(cg_methods);

    colors = {'r', 'b', 'g', 'm'};
    markers = {'o', 's', 'd', '^'};

    % ----- Figure 1: Cost function comparison -----
    figH = figure('Name', 'Total Objective Function', ...
        'Position', [50, 50, 700, 550], 'Color', 'w');
    hold on;
    for i = 1:num_methods
        plot(0:all_results{i}.num_iter-1, all_results{i}.J_history, ...
            [colors{i} '-' markers{i}], 'LineWidth', 1.6, 'MarkerSize', 6, 'MarkerIndices', 1:15:all_results{i}.num_iter);
    end
    hold off;
    set(gca, 'YScale', 'log');   % semilogarithmic (semilogy after hold on stays linear)
    grid on;
    xlabel('iteration');
    ylabel('total cost J');
    title(' ');
    legend(cg_methods, 'Location', 'northeast');
    save_figure_safe(figH, 'cg_comparison_cost_function.png', 600);
    fprintf('Saved cost function plot to cg_comparison_cost_function.png\n');

    % ----- Figure 2: Gradient norm comparison -----
    figH = figure('Name', 'Gradient Norm', ...
        'Position', [50, 50, 700, 550], 'Color', 'w');
    hold on;
    for i = 1:num_methods
        plot(0:all_results{i}.num_iter-1, all_results{i}.grad_norm_history, ...
            [colors{i} '-' markers{i}], 'LineWidth', 1.6, 'MarkerSize', 6, 'MarkerIndices', 1:15:all_results{i}.num_iter);
    end
    hold off;
    set(gca, 'YScale', 'log');   % semilogarithmic (semilogy after hold on stays linear)
    grid on;
    xlabel('iteration');
    ylabel('||\nabla J||');
    title(' ');
    legend(cg_methods, 'Location', 'northeast');
    save_figure_safe(figH, 'cg_comparison_gradient_norm.png', 600);
    fprintf('Saved gradient norm plot to cg_comparison_gradient_norm.png\n');

    % ----- Figure 3: Forward error comparison -----
    figH = figure('Name', 'Forward Error', ...
        'Position', [50, 50, 700, 550], 'Color', 'w');
    hold on;
    for i = 1:num_methods
        plot(0:all_results{i}.num_iter-1, all_results{i}.forward_error_history, ...
            [colors{i} '-' markers{i}], 'LineWidth', 1.6, 'MarkerSize', 6, 'MarkerIndices', 1:15:all_results{i}.num_iter);
    end
    hold off;
    grid on;
    xlabel('iteration');
    ylabel('forward error');
    title(' ');
    legend(cg_methods, 'Location', 'northeast');
    save_figure_safe(figH, 'cg_comparison_forward_error.png', 600);
    fprintf('Saved forward error plot to cg_comparison_forward_error.png\n');

    % ----- Figure 4: Inverse error comparison -----
    figH = figure('Name', 'Inverse Error', ...
        'Position', [50, 50, 700, 550], 'Color', 'w');
    hold on;
    for i = 1:num_methods
        plot(0:all_results{i}.num_iter-1, all_results{i}.inverse_error_history, ...
            [colors{i} '-' markers{i}], 'LineWidth', 1.6, 'MarkerSize', 6, 'MarkerIndices', 1:15:all_results{i}.num_iter);
    end
    hold off;
    set(gca, 'YScale', 'log');   % semilogarithmic (semilogy after hold on stays linear)
    grid on;
    xlabel('iteration');
    ylabel('inverse error');
    title(' ');
    legend(cg_methods, 'Location', 'northeast');
    save_figure_safe(figH, 'cg_comparison_inverse_error.png', 600);
    fprintf('Saved inverse error plot to cg_comparison_inverse_error.png\n');

    % ----- Figure 5: Beta history comparison -----
    figH = figure('Name', 'CG Coefficient Beta_D History', ...
        'Position', [50, 50, 700, 550], 'Color', 'w');
    hold on;
    for i = 1:num_methods
        plot(1:all_results{i}.num_iter, all_results{i}.beta_D_history, ...
            [colors{i} '-' markers{i}], 'LineWidth', 2, 'MarkerSize', 5);
    end
    hold off;
    grid on;
    xlabel('iteration');
    ylabel('\beta_D');
    title('CG Coefficient \beta_D History');
    legend(cg_methods, 'Location', 'best');
    save_figure_safe(figH, 'cg_comparison_beta_history.png', 600);
    fprintf('Saved beta history plot to cg_comparison_beta_history.png\n');

    % ----- Figure 6: CG restart comparison -----
    figH = figure('Name', 'CG Restart Frequency', ...
        'Position', [50, 50, 700, 550], 'Color', 'w');
    restart_counts = zeros(num_methods, 1);
    for i = 1:num_methods
        restart_counts(i) = all_results{i}.cg_restart_count;
    end
    bar(restart_counts);
    set(gca, 'XTickLabel', cg_methods);
    ylabel('number of restarts');
    title('CG Restart Frequency');
    grid on;
    save_figure_safe(figH, 'cg_comparison_restart_frequency.png', 600);
    fprintf('Saved restart frequency plot to cg_comparison_restart_frequency.png\n');
end

%% ============================================================================
%% LINE PROFILES OF THE REFERENCE METHOD
%% ============================================================================
%% Layout per figure:
%%   - 1 row x 3 columns (one panel per species: H2SO4, NaOH, Na2SO4)
%%   - Single color per species: H2SO4=red, NaOH=blue, Na2SO4=green
%%   - 3 line styles for 3 y-stations: y=1 (solid), y=4 (dashed), y=8 (dotted)
%%   - For D plots: horizontal dash-dot reference line (reference values at the
%%     temperature of the data set, from step 1)
%% ============================================================================

function create_line_profile_plots(all_results, cg_methods, X, Y)

    %% ----- Reference method (ref_method in the main script) -----
    best_idx = ref_method_index(all_results, cg_methods);
    best_method   = cg_methods{best_idx};
    R             = all_results{best_idx};

    %% ----- Coordinates -----
    x_vec = X(1, :);    % transverse coordinate x [mm]
    y_vec = Y(:, 1);    % longitudinal coordinate y [mm]

    %% ----- Target y-stations: 1, 4, 8 mm -----
    y_targets = [1, 4, 8];
    n_y       = numel(y_targets);
    y_idx     = zeros(1, n_y);
    y_actual  = zeros(1, n_y);
    for k = 1:n_y
        [~, y_idx(k)] = min(abs(y_vec - y_targets(k)));
        y_actual(k)   = y_vec(y_idx(k));
    end

    %% ----- Styles -----
    species_colors = {[0.85 0.10 0.10], ...   % H2SO4  red
                      [0.10 0.30 0.85], ...   % NaOH   blue
                      [0.10 0.65 0.20]};      % Na2SO4 green
    line_styles    = {'-', '--', ':'};        % y=1, y=4, y=8
    line_width     = 1.6;

    species_names_tex = {'H_2SO_4', 'NaOH', 'Na_2SO_4'};

    unit_D = 'mm^2/s';
    unit_S = 'mol/(L\cdot s)';

    %% ----- Literature D values (mm^2/s) at the temperature of this condition -----
    D_lit = getappdata(0, 'D_lit');

    %% ===================================================================
    %% FIGURE 1 - Source terms S(x) at 3 y-stations
    %% ===================================================================
    fig_S = figure('Name', sprintf('Line Profiles - S(x) [%s]', best_method), ...
                   'Position', [50, 50, 1500, 380], 'Color', 'w');

    S_fields = {'S1_phys', 'S2_phys', 'S3_phys'};

    for sp = 1:3
        subplot(1, 3, sp);
        hold on; grid on; box on;

        F = R.(S_fields{sp});   % Ny x Nx, mol/(L*s)

        h = gobjects(1, n_y);
        for k = 1:n_y
            h(k) = plot(x_vec, F(y_idx(k), :), ...
                'Color',     species_colors{sp}, ...
                'LineStyle', line_styles{k}, ...
                'LineWidth', line_width);
        end
        hold off;

        title(species_names_tex{sp}, 'Interpreter', 'tex', 'FontSize', 9);
        xlabel('x (mm)', 'FontSize', 9);
        ylabel(sprintf('S (%s)', unit_S), 'FontSize', 9, 'Interpreter', 'tex');
        set(gca, 'FontSize', 8);

        % Legend only in first panel
        if sp == 1
            legend(h, ...
                {sprintf('y = %.0f mm', y_actual(1)), ...
                 sprintf('y = %.0f mm', y_actual(2)), ...
                 sprintf('y = %.0f mm', y_actual(3))}, ...
                'Location', 'best', 'FontSize', 7, 'Box', 'off');
        end
    end

    save_figure_safe(fig_S, 'line_profiles_S_transverse.png', 600);
    fprintf('   Saved -> line_profiles_S_transverse.png\n');

    %% ===================================================================
    %% FIGURE 2 - Diffusion coefficients D(x) at 3 y-stations + literature
    %% ===================================================================
    fig_D = figure('Name', sprintf('Line Profiles - D(x) [%s]', best_method), ...
                   'Position', [50, 50, 1500, 380], 'Color', 'w');

    D_fields = {'D1_phys', 'D2_phys', 'D3_phys'};

    for sp = 1:3
        subplot(1, 3, sp);
        hold on; grid on; box on;

        F = R.(D_fields{sp});   % Ny x Nx, mm^2/s

        h = gobjects(1, n_y + 1);
        for k = 1:n_y
            h(k) = plot(x_vec, F(y_idx(k), :), ...
                'Color',     species_colors{sp}, ...
                'LineStyle', line_styles{k}, ...
                'LineWidth', line_width);
        end

        % Literature reference: horizontal dash-dot line
        h(n_y + 1) = yline(D_lit(sp), 'k-.', 'LineWidth', 1.2);

        hold off;

        title(species_names_tex{sp}, 'Interpreter', 'tex', 'FontSize', 9);
        xlabel('x (mm)', 'FontSize', 9);
        ylabel(sprintf('D (%s)', unit_D), 'FontSize', 9);
        set(gca, 'FontSize', 8);

        % Legend only in first panel
        if sp == 1
            legend(h, ...
                {sprintf('y = %.0f mm', y_actual(1)), ...
                 sprintf('y = %.0f mm', y_actual(2)), ...
                 sprintf('y = %.0f mm', y_actual(3)), ...
                 'literature'}, ...
                'Location', 'best', 'FontSize', 7, 'Box', 'off');
        end
    end


    save_figure_safe(fig_D, 'line_profiles_D_transverse.png', 600);
    fprintf('   Saved -> line_profiles_D_transverse.png\n');
end

%% ============================================================================
%% PERFORMANCE COMPARISON TABLE
%% ============================================================================

function create_performance_comparison_table(all_results, cg_methods)
    fprintf('\n');
    fprintf('================================================================================\n');
    fprintf('                          PERFORMANCE COMPARISON TABLE                          \n');
    fprintf('================================================================================\n');
    fprintf('\n');

    % Table header
    fprintf('%-10s | %12s | %10s | %12s | %12s | %12s | %10s\n', ...
        'Method', 'Final Cost', 'Iterations', 'Grad Norm', 'Forward Err', 'Inverse Err', 'Time (s)');
    fprintf('%s\n', repmat('-', 1, 100));

    % Table data
    for i = 1:length(cg_methods)
        fprintf('%-10s | %.6e | %10d | %.6e | %.6e | %.6e | %10.2f\n', ...
            cg_methods{i}, ...
            all_results{i}.best_cost, ...
            all_results{i}.num_iter, ...
            all_results{i}.final_grad_norm, ...
            all_results{i}.forward_error_history(end), ...
            all_results{i}.inverse_error_history(end), ...
            all_results{i}.total_time);
    end

    fprintf('%s\n', repmat('-', 1, 100));
    fprintf('\n');

    % CG Method formulas
    fprintf('CG Method Formulas:\n');
    fprintf('  FR (Fletcher-Reeves):  beta = ||g_k||^2 / ||g_{k-1}||^2\n');
    fprintf('  PR (Polak-Ribiere+):   beta = max(0, g_k*(g_k - g_{k-1}) / ||g_{k-1}||^2)\n');
    fprintf('  PB (Powell-Beale):     beta = g_k*(g_k - g_{k-1}) / (d_{k-1}*(g_k - g_{k-1}))\n');
    fprintf('                         with restart when |g_k*g_{k-1}| >= 0.2||g_k||^2\n');
    fprintf('\n');

    % Additional statistics
    fprintf('Additional Statistics:\n');
    fprintf('%-10s | %12s | %15s | %18s\n', 'Method', 'CG Restarts', 'Best Iter', 'Avg Time/Iter (s)');
    fprintf('%s\n', repmat('-', 1, 70));

    for i = 1:length(cg_methods)
        avg_time_per_iter = all_results{i}.total_time / all_results{i}.num_iter;
        fprintf('%-10s | %12d | %15d | %18.3f\n', ...
            cg_methods{i}, ...
            all_results{i}.cg_restart_count, ...
            all_results{i}.best_iter, ...
            avg_time_per_iter);
    end

    fprintf('%s\n', repmat('-', 1, 70));
    fprintf('\n');
end

%% ============================================================================
%% PHYSICAL PARAMETER COMPARISON PLOTS (3D SURFACES)
%% ============================================================================

function create_parameter_comparison_plots_physical(all_results, cg_methods, X, Y)
    num_methods = length(cg_methods);

    % --- Sanitize any non-finite values that can crash surf/export ---
    flds = {'D1_phys','D2_phys','D3_phys','S1_phys','S2_phys','S3_phys'};
    for i = 1:num_methods
        for f = 1:numel(flds)
            v = all_results{i}.(flds{f});
            v(~isfinite(v)) = 0;
            all_results{i}.(flds{f}) = v;
        end
    end

    LX = 'x (mm)';
    LY = 'y (mm)';
    unit_D = 'mm^2/s';
    unit_S = 'mol/(L\cdot s)';

    % Physical Diffusion coefficient D1 comparison
    figH = figure('Name', 'Physical Parameter Comparison - D1_phys', ...
           'Position', [50, 50, 1400, 400], 'Color', 'w');
    for i = 1:num_methods
        subplot(1, num_methods, i);
        D1_phys = all_results{i}.D1_phys;
        surf(X, Y, D1_phys, 'EdgeColor', 'none');
        colorbar;
        colormap('parula');
        title(sprintf('D_1 - %s', cg_methods{i}), 'Interpreter', 'tex');
        xlabel(LX); ylabel(LY); zlabel(sprintf('D_1 (%s)', unit_D));
        view(45, 30);
    end
    save_figure_safe(figH, 'parameter_comparison_D1_phys.png', 600);

    % Physical Diffusion coefficient D2 comparison
    figH = figure('Name', 'Physical Parameter Comparison - D2_phys', ...
           'Position', [50, 50, 1400, 400], 'Color', 'w');
    for i = 1:num_methods
        subplot(1, num_methods, i);
        D2_phys = all_results{i}.D2_phys;
        surf(X, Y, D2_phys, 'EdgeColor', 'none');
        colorbar;
        colormap('parula');
        title(sprintf('D_2 - %s', cg_methods{i}), 'Interpreter', 'tex');
        xlabel(LX); ylabel(LY); zlabel(sprintf('D_2 (%s)', unit_D));
        view(45, 30);
    end
    save_figure_safe(figH, 'parameter_comparison_D2_phys.png', 600);

    % Physical Diffusion coefficient D3 comparison
    figH = figure('Name', 'Physical Parameter Comparison - D3_phys', ...
           'Position', [50, 50, 1400, 400], 'Color', 'w');
    for i = 1:num_methods
        subplot(1, num_methods, i);
        D3_phys = all_results{i}.D3_phys;
        surf(X, Y, D3_phys, 'EdgeColor', 'none');
        colorbar;
        colormap('parula');
        title(sprintf('D_3 - %s', cg_methods{i}), 'Interpreter', 'tex');
        xlabel(LX); ylabel(LY); zlabel(sprintf('D_3 (%s)', unit_D));
        view(45, 30);
    end
    save_figure_safe(figH, 'parameter_comparison_D3_phys.png', 600);

    % Physical Source term S1 comparison
    figH = figure('Name', 'Physical Parameter Comparison - S1_phys', ...
           'Position', [50, 50, 1400, 400], 'Color', 'w');
    for i = 1:num_methods
        subplot(1, num_methods, i);
        S1_phys = all_results{i}.S1_phys;
        surf(X, Y, S1_phys, 'EdgeColor', 'none');
        colorbar;
        colormap('parula');
        title(sprintf('S_1 - %s', cg_methods{i}), 'Interpreter', 'tex');
        xlabel(LX); ylabel(LY); zlabel(sprintf('S_1 (%s)', unit_S));
        view(45, 30);
    end
    save_figure_safe(figH, 'parameter_comparison_S1_phys.png', 600);

    % Physical Source term S2 comparison
    figH = figure('Name', 'Physical Parameter Comparison - S2_phys', ...
           'Position', [50, 50, 1400, 400], 'Color', 'w');
    for i = 1:num_methods
        subplot(1, num_methods, i);
        S2_phys = all_results{i}.S2_phys;
        surf(X, Y, S2_phys, 'EdgeColor', 'none');
        colorbar;
        colormap('parula');
        title(sprintf('S_2 - %s', cg_methods{i}), 'Interpreter', 'tex');
        xlabel(LX); ylabel(LY); zlabel(sprintf('S_2 (%s)', unit_S));
        view(45, 30);
    end
    save_figure_safe(figH, 'parameter_comparison_S2_phys.png', 600);

    % Physical Source term S3 comparison
    figH = figure('Name', 'Physical Parameter Comparison - S3_phys', ...
           'Position', [50, 50, 1400, 400], 'Color', 'w');
    for i = 1:num_methods
        subplot(1, num_methods, i);
        S3_phys = all_results{i}.S3_phys;
        surf(X, Y, S3_phys, 'EdgeColor', 'none');
        colorbar;
        colormap('parula');
        title(sprintf('S_3 - %s', cg_methods{i}), 'Interpreter', 'tex');
        xlabel(LX); ylabel(LY); zlabel(sprintf('S_3 (%s)', unit_S));
        view(45, 30);
    end
    save_figure_safe(figH, 'parameter_comparison_S3_phys.png', 600);

    fprintf('Saved physical parameter comparison plots\n');

    %% BEST METHOD SUMMARY: 3D shared-axis figures (D and S)
    best_idx = ref_method_index(all_results, cg_methods);
    best_cost_val = all_results{best_idx}.best_cost;
    best_method = cg_methods{best_idx};

    D1p = all_results{best_idx}.D1_phys;
    D2p = all_results{best_idx}.D2_phys;
    D3p = all_results{best_idx}.D3_phys;
    S1p = all_results{best_idx}.S1_phys;
    S2p = all_results{best_idx}.S2_phys;
    S3p = all_results{best_idx}.S3_phys;

    %% Shared scales
    D_all   = [D1p(:); D2p(:); D3p(:)];
    D_min   = min(D_all);   D_max = max(D_all);

    S_all   = [S1p(:); S2p(:); S3p(:)];
    S_min   = min(S_all);   S_max = max(S_all);

    %% FIGURE: Diffusion Coefficients 3D
    figH = figure('Name', sprintf('Physical Diffusion Coefficients - Reference Method: %s', best_method), ...
    'Position', [100, 100, 1400, 500], 'Color', 'w');

    subplot(1, 3, 1);
    surf(X, Y, D1p, 'EdgeColor', 'none', 'FaceColor', 'interp');
    colormap(gca, 'parula');
    caxis([D_min D_max]); zlim([D_min D_max]);
    xlim([min(X(:)) max(X(:))]); ylim([min(Y(:)) max(Y(:))]);
    title(sprintf('D_1 (%s)', best_method), 'Interpreter', 'tex');
    xlabel(LX); ylabel(LY); zlabel(sprintf('D_1 (%s)', unit_D));
    view(45, 30); lighting gouraud; camlight('headlight'); grid on;

    subplot(1, 3, 2);
    surf(X, Y, D2p, 'EdgeColor', 'none', 'FaceColor', 'interp');
    colormap(gca, 'parula');
    caxis([D_min D_max]); zlim([D_min D_max]);
    xlim([min(X(:)) max(X(:))]); ylim([min(Y(:)) max(Y(:))]);
    title(sprintf('D_2 (%s)', best_method), 'Interpreter', 'tex');
    xlabel(LX); ylabel(LY); zlabel(sprintf('D_2 (%s)', unit_D));
    view(45, 30); lighting gouraud; camlight('headlight'); grid on;

    subplot(1, 3, 3);
    surf(X, Y, D3p, 'EdgeColor', 'none', 'FaceColor', 'interp');
    colormap(gca, 'parula');
    caxis([D_min D_max]); zlim([D_min D_max]);
    xlim([min(X(:)) max(X(:))]); ylim([min(Y(:)) max(Y(:))]);
    title(sprintf('D_3 (%s)', best_method), 'Interpreter', 'tex');
    xlabel(LX); ylabel(LY); zlabel(sprintf('D_3 (%s)', unit_D));
    view(45, 30); lighting gouraud; camlight('headlight'); grid on;

    c = colorbar('Position', [0.93 0.15 0.015 0.7]);
    ylabel(c, sprintf('D (%s)', unit_D), 'FontSize', 10);

    save_figure_safe(figH, 'physical_parameters_summary_D_best.png', 600);
    fprintf('Saved 3D diffusion summary plot (reference method: %s)\n', best_method);

    %% FIGURE: Source Terms 3D
    figH = figure('Name', sprintf('Physical Source Terms - Reference Method: %s', best_method), ...
    'Position', [100, 100, 1400, 500], 'Color', 'w');

    subplot(1, 3, 1);
    surf(X, Y, S1p, 'EdgeColor', 'none', 'FaceColor', 'interp');
    colormap(gca, 'parula');
    caxis([S_min S_max]); zlim([S_min S_max]);
    xlim([min(X(:)) max(X(:))]); ylim([min(Y(:)) max(Y(:))]);
    title(sprintf('S_1 (%s)', best_method), 'Interpreter', 'tex');
    xlabel(LX); ylabel(LY); zlabel(sprintf('S_1 (%s)', unit_S));
    view(45, 30); lighting gouraud; camlight('headlight'); grid on;

    subplot(1, 3, 2);
    surf(X, Y, S2p, 'EdgeColor', 'none', 'FaceColor', 'interp');
    colormap(gca, 'parula');
    caxis([S_min S_max]); zlim([S_min S_max]);
    xlim([min(X(:)) max(X(:))]); ylim([min(Y(:)) max(Y(:))]);
    title(sprintf('S_2 (%s)', best_method), 'Interpreter', 'tex');
    xlabel(LX); ylabel(LY); zlabel(sprintf('S_2 (%s)', unit_S));
    view(45, 30); lighting gouraud; camlight('headlight'); grid on;

    subplot(1, 3, 3);
    surf(X, Y, S3p, 'EdgeColor', 'none', 'FaceColor', 'interp');
    colormap(gca, 'parula');
    caxis([S_min S_max]); zlim([S_min S_max]);
    xlim([min(X(:)) max(X(:))]); ylim([min(Y(:)) max(Y(:))]);
    title(sprintf('S_3 (%s)', best_method), 'Interpreter', 'tex');
    xlabel(LX); ylabel(LY); zlabel(sprintf('S_3 (%s)', unit_S));
    view(45, 30); lighting gouraud; camlight('headlight'); grid on;

    c = colorbar('Position', [0.93 0.15 0.015 0.7]);
    ylabel(c, sprintf('S (%s)', unit_S), 'FontSize', 10);

    save_figure_safe(figH, 'physical_parameters_summary_S_best.png', 600);
    fprintf('Saved 3D source terms summary plot (reference method: %s)\n', best_method);

    %% FIGURE: 2x3 Combined Summary
    figH = figure('Name', 'Physical Parameters Summary - Reference Method', ...
        'Position', [50, 50, 1650, 900], 'Color', 'w');

    %% Shared view angle for all 3D subplots
    az = 35;
    el = 28;

    param_data   = {D1p,  D2p,  D3p,  S1p,  S2p,  S3p};
    param_titles = { ...
        sprintf('D_1 (%s)', best_method), ...
        sprintf('D_2 (%s)', best_method), ...
        sprintf('D_3 (%s)', best_method), ...
        sprintf('S_1 (%s)', best_method), ...
        sprintf('S_2 (%s)', best_method), ...
        sprintf('S_3 (%s)', best_method)};
    param_zlabels = { ...
        sprintf('D_1 (%s)', unit_D), sprintf('D_2 (%s)', unit_D), sprintf('D_3 (%s)', unit_D), ...
        sprintf('S_1 (%s)', unit_S), sprintf('S_2 (%s)', unit_S), sprintf('S_3 (%s)', unit_S)};
    clims = {[D_min D_max],[D_min D_max],[D_min D_max], ...
             [S_min S_max],[S_min S_max],[S_min S_max]};
    cmaps = {'parula','parula','parula','hot','hot','hot'};

    for k = 1:6
        ax = subplot(2, 3, k);
        Fp = param_data{k};

        surf(X, Y, Fp, 'EdgeColor', 'none', 'FaceAlpha', 0.92);
        colormap(ax, cmaps{k});
        if clims{k}(2) > clims{k}(1)
            caxis(clims{k});
        end

        cb = colorbar('Location', 'eastoutside');
        cb.Label.String  = param_zlabels{k};
        cb.Label.FontSize = 9;

        xlabel(LX, 'FontSize', 9);
        ylabel(LY, 'FontSize', 9);
        zlabel(param_zlabels{k}, 'FontSize', 9);
        title(param_titles{k}, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'tex');

        % Min/max annotation inset
        text(0.02, 0.97, sprintf('min = %.3e\nmax = %.3e', min(Fp(:)), max(Fp(:))), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'FontSize', 7.5, 'Color', [0.15 0.15 0.15], ...
            'BackgroundColor', [1 1 1 0.65], 'EdgeColor', [0.7 0.7 0.7]);

        view(az, el);
        grid on; box on;
        set(ax, 'FontSize', 8);
        lighting phong;
        camlight('headlight');
    end

    sgtitle(sprintf('Physical Parameters - Reference Method: %s   (J = %.4e)', ...
        best_method, best_cost_val), ...
        'FontSize', 13, 'FontWeight', 'bold');

    save_figure_safe(figH, 'physical_parameters_summary_best.png', 600);
    fprintf('Saved physical parameters summary plot (3D)\n');
end

%% ============================================================================
%% SAVE COMPARISON RESULTS
%% ============================================================================

function save_comparison_results(all_results, cg_methods, X, Y, x_phys, y_phys, ...
    alpha_D_reg, alpha_S_reg, reg_type, param_bounds, D_phys, S_phys, ndD, ndS)

    comparison = struct();
    comparison.cg_methods = cg_methods;
    comparison.all_results = all_results;
    comparison.X = X;
    comparison.Y = Y;
    comparison.x_phys = x_phys;
    comparison.y_phys = y_phys;
    comparison.alpha_D_reg = alpha_D_reg;
    comparison.alpha_S_reg = alpha_S_reg;
    comparison.reg_type = reg_type;
    comparison.param_bounds = param_bounds;
    comparison.comparison_date = datestr(now);
    comparison.D_phys = D_phys;
    comparison.S_phys = S_phys;
    comparison.ndD = ndD;
    comparison.ndS = ndS;
    comparison.implementation_note = 'TRUE Powell-Beale: beta = g_k*(g_k-g_{k-1}) / (d_{k-1}*(g_k-g_{k-1}))';

    save('cg_method_comparison_results.mat', 'comparison');
    fprintf('Saved comparison results to cg_method_comparison_results.mat\n');

    save('physical_parameters_results.mat', 'D_phys', 'S_phys', 'ndD', 'ndS', 'X', 'Y', 'cg_methods');
    fprintf('Saved physical parameters to physical_parameters_results.mat\n');

    % Write detailed comparison report
    fid = fopen('cg_method_comparison_report.txt', 'w');
    fprintf(fid, '===============================================================\n');
    fprintf(fid, '  CONJUGATE GRADIENT METHOD COMPARISON REPORT\n');
    fprintf(fid, '  (with TRUE Powell-Beale Implementation)\n');
    fprintf(fid, '===============================================================\n\n');
    fprintf(fid, 'Date: %s\n\n', datestr(now));

    fprintf(fid, 'Methods Compared:\n');
    fprintf(fid, '  1. FR - Fletcher-Reeves\n');
    fprintf(fid, '     beta = ||g_k||^2 / ||g_{k-1}||^2\n\n');
    fprintf(fid, '  2. PR - Polak-Ribiere (PR+)\n');
    fprintf(fid, '     beta = max(0, g_k*(g_k - g_{k-1}) / ||g_{k-1}||^2)\n\n');
    fprintf(fid, '  3. PB - Powell-Beale (TRUE)\n');
    fprintf(fid, '     beta = g_k*(g_k - g_{k-1}) / (d_{k-1}*(g_k - g_{k-1}))\n');
    fprintf(fid, '     with restart when |g_k*g_{k-1}| >= 0.2||g_k||^2\n\n');

    fprintf(fid, 'Non-dimensionalization Parameters:\n');
    fprintf(fid, '  ndD (Diffusion scaling): %.4e mm^2/s\n', ndD);
    fprintf(fid, '  ndS (Source scaling):    %.4e mol/L/s\n\n', ndS);

    fprintf(fid, 'PERFORMANCE SUMMARY\n');
    fprintf(fid, '%s\n', repmat('-', 1, 60));

    for i = 1:length(cg_methods)
        fprintf(fid, '\n%s Method:\n', cg_methods{i});
        fprintf(fid, '  Final Cost:        %.6e\n', all_results{i}.best_cost);
        fprintf(fid, '  Iterations:        %d\n', all_results{i}.num_iter);
        fprintf(fid, '  Gradient Norm:     %.6e\n', all_results{i}.final_grad_norm);
        fprintf(fid, '  CG Restarts:       %d\n', all_results{i}.cg_restart_count);
        fprintf(fid, '  Computation Time:  %.2f seconds\n', all_results{i}.total_time);

        fprintf(fid, '\n  Physical Parameter Ranges:\n');
        fprintf(fid, '    D1_phys: [%.4e, %.4e] mm^2/s\n', ...
            min(all_results{i}.D1_phys(:)), max(all_results{i}.D1_phys(:)));
        fprintf(fid, '    D2_phys: [%.4e, %.4e] mm^2/s\n', ...
            min(all_results{i}.D2_phys(:)), max(all_results{i}.D2_phys(:)));
        fprintf(fid, '    D3_phys: [%.4e, %.4e] mm^2/s\n', ...
            min(all_results{i}.D3_phys(:)), max(all_results{i}.D3_phys(:)));
        fprintf(fid, '    S1_phys: [%.4e, %.4e] mol/L/s\n', ...
            min(all_results{i}.S1_phys(:)), max(all_results{i}.S1_phys(:)));
        fprintf(fid, '    S2_phys: [%.4e, %.4e] mol/L/s\n', ...
            min(all_results{i}.S2_phys(:)), max(all_results{i}.S2_phys(:)));
        fprintf(fid, '    S3_phys: [%.4e, %.4e] mol/L/s\n', ...
            min(all_results{i}.S3_phys(:)), max(all_results{i}.S3_phys(:)));
    end

    fprintf(fid, '\n===============================================================\n');
    fclose(fid);

    fprintf('Saved detailed comparison report to cg_method_comparison_report.txt\n');
end

%% ============================================================================
%% DISPLAY METHOD RECOMMENDATION
%% ============================================================================

function display_method_recommendation(all_results, cg_methods)
    fprintf('\n');
    fprintf('================================================================\n');
    fprintf('                    RECOMMENDATION SUMMARY                      \n');
    fprintf('              (TRUE Powell-Beale Implementation)                \n');
    fprintf('================================================================\n');
    fprintf('\n');

    % Find best by different criteria
    [best_cost, best_cost_idx] = min(cellfun(@(x) x.best_cost, all_results));
    [best_time, best_time_idx] = min(cellfun(@(x) x.total_time, all_results));
    [best_iter, best_iter_idx] = min(cellfun(@(x) x.num_iter, all_results));
    [min_restarts, min_restart_idx] = min(cellfun(@(x) x.cg_restart_count, all_results));

    fprintf('BEST BY CRITERION:\n\n');
    fprintf('  * Lowest objective J: %s (J: %.6e)\n', ...
        cg_methods{best_cost_idx}, best_cost);
    fprintf('  * Computation Speed: %s (Time: %.2f s)\n', ...
        cg_methods{best_time_idx}, best_time);
    fprintf('  * Efficiency:        %s (Iterations: %d)\n', ...
        cg_methods{best_iter_idx}, best_iter);
    fprintf('  * Stability:         %s (CG Restarts: %d)\n', ...
        cg_methods{min_restart_idx}, min_restarts);

    ri = ref_method_index(all_results, cg_methods);
    fprintf('\n  Reference method used for the fields and figures: %s (set by ref_method).\n', cg_methods{ri});
    fprintf('  A lower J does not imply more accurate fields: compare the data misfit and the synthetic test.\n');
    fprintf('\n');
    fprintf('CG METHOD FORMULAS:\n\n');
    fprintf('  FR: beta = ||g_k||^2 / ||g_{k-1}||^2\n');
    fprintf('  PR: beta = max(0, g_k*(g_k - g_{k-1}) / ||g_{k-1}||^2)\n');
    fprintf('  PB: beta = g_k*(g_k - g_{k-1}) / (d_{k-1}*(g_k - g_{k-1}))\n');
    fprintf('      with restart when |g_k*g_{k-1}| >= 0.2||g_k||^2\n');
    fprintf('  DY: beta = ||g_k||^2 / (d_{k-1}*(g_k - g_{k-1}))\n');
    fprintf('\n');

    % Overall recommendation
    if best_cost_idx == best_time_idx && best_cost_idx == best_iter_idx
        fprintf('OVERALL RECOMMENDATION: %s\n', cg_methods{best_cost_idx});
        fprintf('   This method dominates all criteria.\n');
    else
        fprintf('RECOMMENDATION DEPENDS ON PRIORITY:\n');
        fprintf('   * For the lowest objective J -> Use %s\n', cg_methods{best_cost_idx});
        fprintf('   * For fastest results -> Use %s\n', cg_methods{best_time_idx});
        fprintf('   * For most efficient -> Use %s\n', cg_methods{best_iter_idx});
        fprintf('   * For most stable -> Use %s\n', cg_methods{min_restart_idx});
    end

    fprintf('\n');
    fprintf('PHYSICAL PARAMETERS AVAILABLE IN WORKSPACE:\n');
    fprintf('   * D_phys.FR, D_phys.PR, D_phys.PB (structures with D1, D2, D3)\n');
    fprintf('   * S_phys.FR, S_phys.PR, S_phys.PB (structures with S1, S2, S3)\n');
    fprintf('   * D1_phys_best, D2_phys_best, D3_phys_best (from reference method)\n');
    fprintf('   * S1_phys_best, S2_phys_best, S3_phys_best (from reference method)\n');
    fprintf('\n');
    fprintf('OUTPUT FIGURES GENERATED:\n');
    fprintf('   * line_profiles_D_transverse.png           (D(x) reference format)\n');
    fprintf('   * line_profiles_S_transverse.png           (S(x) reference format)\n');
    fprintf('   * parameter_comparison_D1/D2/D3_phys.png   (3D surfaces, per-method)\n');
    fprintf('   * parameter_comparison_S1/S2/S3_phys.png   (3D surfaces, per-method)\n');
    fprintf('   * physical_parameters_summary_D_best.png   (3D, best D)\n');
    fprintf('   * physical_parameters_summary_S_best.png   (3D, best S)\n');
    fprintf('   * physical_parameters_summary_best.png     (2x3 summary)\n');
    fprintf('\n');
    fprintf('================================================================\n');
    fprintf('\n');
end
