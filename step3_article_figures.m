%STEP3_ARTICLE_FIGURES  All data figures and numbers of the article from the saved results.
%
%   Step 3 of 3 (optional). Reads the results of steps 1 and 2 and draws every data
%   figure of the article (Figs. 2(b)-13) in the layout and scales of the article.
%   No inversion is repeated; runs in well under a minute.
%
%   Input    results/<condition>/step1_initial_estimates.mat  (from step 1)
%            results/<condition>/step2_cg_results.mat         (from step 2)
%            results/<condition>/map_concentration_medfilt.mat (from step 1, Fig. 2(b))
%            results/<condition>/step2_lcurve_results.mat     (optional, Fig. 5; see below)
%            results/T20, T28, T37/...                         (Figs. 12-13; see below)
%   Output   in results/<condition>/:
%            Fig02b_measured_maps.png    measured maps (3 x 3 median filter)      Fig. 2(b)
%            Fig02c_fitted_maps.png      fitted concentration fields              Fig. 2(c)
%            Fig03_fit_quality.png       measured vs fitted profiles             Fig. 3
%            Fig04_direct_estimates.png  initial D (direct) and S (balance)      Fig. 4
%            Fig05_L_curves.png          L-curves of the regularization weights  Fig. 5
%            Fig06_convergence.png       (a) total cost J, (b) gradient norm     Fig. 6
%            Fig07_errors.png            (a) forward error, (b) inverse error    Fig. 7
%            Fig08a_D_3D_<m>.png         D fields, 3D, common scale              Fig. 8(a)
%            Fig08b_D_lines_<m>.png      D profiles at y = 1, 4, 8 mm            Fig. 8(b)
%            Fig09_seq_S_<m>.png         simultaneous vs sequential recovery, S  Fig. 9
%            Fig10a_S_3D_<m>.png         S fields, 3D, common scale              Fig. 10(a)
%            Fig10b_S_lines_<m>.png      S profiles at y = 1, 4, 8 mm            Fig. 10(b)
%            Fig11_MC_std.png            Monte Carlo standard deviations         Fig. 11
%            step3_article_numbers.txt   numbers for Table 2 and Sections 6.2-6.4
%            (<m> = CG method, FR by default)
%            in results/:
%            Fig12_D_temperature.png     D at 20, 28, 37 degC, y = 4 mm          Fig. 12
%            Fig13_S_temperature.png     S at 20, 28, 37 degC, y = 4 mm          Fig. 13
%
%   Fig. 1 and Fig. 2(a) are drawings of the set-up and are not produced by the code.
%
%   Fig. 5   The L-curves need about 15 extra inversions. They are taken, in this order, from
%            (1) results/<condition>/step2_lcurve_results.mat, written by step 2 with
%                alpha_mode = 'lcurve2' and lcurve_only = true;
%            (2) LC/LC2 in step2_cg_results.mat (step 2 run with alpha_mode 'lcurve2');
%            (3) data/<condition>/lcurve_article.csv, the L-curve values of the article
%                run (300 iterations, provided for T28 only).
%   Figs. 12-13  need steps 1 and 2 for T20, T28 and T37 (300 iterations for all
%            data sets, as in the article). The data sets that have been run are
%            drawn; missing ones are listed in the command window.
%
%   Usage    Run after steps 1 and 2 with the same condition.
%   Requires MATLAB R2020a or newer (base MATLAB only)
%
%   Code accompanying the article
%     R. Abdurrahman, K. Matsunaga, Y. Abe, N. Kakuta, "Simultaneous Regularized
%     Field Recovery of Diffusion Coefficients and Source Terms in Microfluidic
%     Acid-Base Reactions",
%     Chemical Engineering Science (submitted; DOI to be added on publication).
%   Equation, figure and table numbers in the comments refer to this article.
%
%   Repository : https://github.com/ridwan029-dev/ychannel-adjoint-inverse (v1.0)
%   Author     : Ridwan Abdurrahman, Kakuta Laboratory, Department of Mechanical
%                Systems Engineering, Tokyo Metropolitan University
%   License    : MIT (see LICENSE)

clear; clc; close all;

condition  = 'T28';                      % data set, as in steps 1 and 2
fig_method = 'FR';                       % method shown in the article (reference method of step 2)
y_lines    = [1 4 8];                    % mm, stations of the line profiles
species    = {'H_2SO_4', 'NaOH', 'Na_2SO_4'};
names      = {'H2SO4', 'NaOH', 'Na2SO4'};
cols       = {[0.82 0.11 0.11], [0.12 0.31 0.85], [0.08 0.60 0.21]};   % species colours
ycols      = {[0.00 0.45 0.70], [0.84 0.37 0.00], [0.00 0.62 0.45]};   % y = 1, 4, 8 mm (Figs. 3-4)
lstyle     = {'-', '--', ':'};           % y = 1, 4, 8 mm
ymark      = {'o', 's', '^'};
clim_C     = [1.2 2.4 0.6];              % mol/L, colour limits of the concentration maps (Fig. 2)
scD        = 1e3;                        % D plotted in 10^-3 mm^2/s
scS        = 1e2;                        % S plotted in 10^-2 mol/(L s)
x_int      = 0.05;                       % mm, |x| < x_int: interfacial median of D
x_max_int  = 0.1;                        % mm, |x| < x_max_int: interfacial maximum of D
x_bulk     = 0.25;                       % mm, outside the reaction zone: x < -x_bulk (H2SO4),
                                         % x > x_bulk (NaOH), |x| > x_bulk (Na2SO4) (Table 2)
t_conds    = {'T20', 'T28', 'T37'};      % Figs. 12-13: temperature series at 0.15 mm/s
t_cols     = {[0.12 0.37 0.75], [0.16 0.62 0.29], [0.75 0.22 0.17]};
t_ls       = {'-', '--', '-.'};
t_mark     = {'o', 's', '^'};
t_y        = 4;                          % mm, station of Figs. 12-13
fs         = 11;                         % font size

start_dir = pwd;
repo_root = fileparts(mfilename('fullpath'));
if isempty(repo_root), repo_root = start_dir; end
res_dir = fullfile(repo_root, 'results', condition);
assert(exist(fullfile(res_dir, 'step2_cg_results.mat'), 'file') == 2, ...
       'Run steps 1 and 2 with condition = %s first (no results in %s).', condition, res_dir);
cd(res_dir);
L = load('step1_initial_estimates.mat');
D_lit = L.D_lit_array;                   % mm^2/s, reference values at the temperature of this data set

R = load('step2_cg_results.mat');
cg_methods = cellstr(R.cg_methods);
m = find(strcmpi(cg_methods, fig_method), 1);
assert(~isempty(m), 'Method %s not found in step2_cg_results.mat', fig_method);
res = R.all_results{m};
x = R.x(:).';  y = R.y(:);  X = R.X;  Y = R.Y;
D = {res.D1_phys, res.D2_phys, res.D3_phys};
S = {res.S1_phys, res.S2_phys, res.S3_phys};
iy = arrayfun(@(v) nearest_index(y, v), y_lines);
ylab = arrayfun(@(v) sprintf('y = %d mm', v), y_lines, 'UniformOutput', false);
done = {};  skipped = {};

fid = fopen('step3_article_numbers.txt', 'w');
pr = @(varargin) [fprintf(varargin{:}), fprintf(fid, varargin{:})];
pr('Numbers of the article, condition %s, method %s (%s)\n\n', condition, fig_method, datestr(now));

%% Fig. 2(b),(c): measured (median-filtered) and fitted concentration maps
Cfit = {};
if isfield(R, 'C_fit'), Cfit = R.C_fit; elseif isfield(L, 'C1'), Cfit = {L.C1, L.C2, L.C3}; end
if exist('map_concentration_medfilt.mat', 'file') == 2
    Mf = load('map_concentration_medfilt.mat');
    Cm_map = {Mf.C1_med, Mf.C2_med, Mf.C3_med};  xm = Mf.x_phys_mm(:).';  ym = Mf.y_phys_mm(:).';
else
    Cm_map = cellfun(@medfilt3x3, {L.C1_obs, L.C2_obs, L.C3_obs}, 'UniformOutput', false);  xm = x;  ym = y.';
end
map_sets = {Cm_map, Cfit};  map_x = {xm, x};  map_y = {ym, y.'};
map_file = {'Fig02b_measured_maps.png', 'Fig02c_fitted_maps.png'};
for q = 1:2
    if isempty(map_sets{q}), continue; end
    figH = figure('Color', 'w', 'Position', [60 60 1500 400]);
    for s = 1:3
        ax = axes('Position', [0.05 + (s - 1)*0.325, 0.42, 0.27, 0.48]);
        imagesc(map_y{q}, map_x{q}, map_sets{q}{s}.');  axis xy;  hold on;
        if q == 1
            for k = 1:3, plot(y_lines(k)*[1 1], [-0.5 0.5], 'w--', 'LineWidth', 1.2); end
        end
        hold off;  caxis([0 clim_C(s)]);  colormap(ax, parula);
        ylim([-0.5 0.5]);  set(ax, 'YTick', [-0.5 0 0.5], 'FontSize', fs, 'Layer', 'top');
        xlabel('y (mm)');  title(species{s}, 'Interpreter', 'tex', 'FontWeight', 'normal');
        if s == 1, ylabel('x (mm)'); end
        cb = colorbar(ax, 'southoutside');
        set(cb, 'Position', [0.05 + (s - 1)*0.325, 0.16, 0.27, 0.045]);
        cbar_ticks(cb, [0 0.5 1]*clim_C(s));  cbar_label(cb, 'C (mol/L)');
    end
    save_figure_safe(figH, map_file{q}, 300);  done{end+1} = map_file{q}; %#ok<SAGROW>
end

%% Fig. 3: measured vs fitted transverse profiles at y = 1, 4, 8 mm
if ~isempty(Cfit)
    Cm = {L.C1_obs, L.C2_obs, L.C3_obs};
    figH = figure('Color', 'w', 'Position', [60 60 1500 450]);
    for s = 1:3
        subplot(1, 3, s); hold on;
        h = gobjects(1, 3);
        for k = 1:3
            plot(x, Cm{s}(iy(k), :), ymark{k}, 'Color', ycols{k}, 'MarkerSize', 5, 'LineWidth', 1);
            plot(x, Cfit{s}(iy(k), :), '-', 'Color', ycols{k}, 'LineWidth', 2);
            h(k) = plot(NaN, NaN, ['-' ymark{k}], 'Color', ycols{k}, 'LineWidth', 2, 'MarkerSize', 5);
        end
        hold off; grid on; box on; xlim([-0.5 0.5]); set(gca, 'XTick', -0.5:0.25:0.5, 'FontSize', fs);
        xlabel('x (mm)'); title(species{s}, 'Interpreter', 'tex', 'FontWeight', 'normal');
        if s == 1
            ylabel('C (mol/L)');
            legend(h, ylab, 'Location', 'northeast');
        end
    end
    save_figure_safe(figH, 'Fig03_fit_quality.png', 300);  done{end+1} = 'Fig03_fit_quality.png';
end

%% Fig. 4: initial estimates of step 1, (a) D from Eq. (20), (b) S from the integral balance
if isfield(L, 'D1_refined')
    D0 = {L.D1_refined, L.D2_refined, L.D3_refined};  S0 = {L.S1, L.S2, L.S3};
    figH = figure('Color', 'w', 'Position', [60 60 1000 1000]);
    for s = 1:3
        for q = 1:2
            subplot(3, 2, 2*(s - 1) + q); hold on;
            h = gobjects(1, 4);
            for k = 1:3
                if q == 1, f = D0{s}(iy(k), :)*scD; else, f = S0{s}(iy(k), :)*scS; end
                h(k) = plot_lm(x, f, ycols{k}, lstyle{k}, ymark{k}, 5);
            end
            if q == 1
                h(4) = plot([-0.5 0.5], D_lit(s)*scD*[1 1], 'k-.', 'LineWidth', 1.1);
                ylabel(sprintf('D of %s (10^{-3} mm^2/s)', species{s}), 'Interpreter', 'tex');
                if s == 1, legend(h, [ylab, {'reference'}], 'Location', 'northwest'); end
            else
                plot([-0.5 0.5], [0 0], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
                ylabel(sprintf('S of %s (10^{-2} mol/(L s))', species{s}), 'Interpreter', 'tex');
            end
            hold off; grid on; box on; xlim([-0.5 0.5]); set(gca, 'XTick', -0.5:0.25:0.5, 'FontSize', fs);
            if s == 1, title(char('(a)' + (q - 1)*[0 1 0]), 'FontWeight', 'normal'); end
            if s == 3, xlabel('x (mm)'); end
        end
    end
    save_figure_safe(figH, 'Fig04_direct_estimates.png', 300);  done{end+1} = 'Fig04_direct_estimates.png';
else
    skipped{end+1} = 'Fig. 4 (no D1_refined in step1_initial_estimates.mat; rerun step 1)';
end

%% Fig. 5: L-curves  (a) common weight, (b) alphaD sweep, (c) alphaS sweep
[T5, sel5, src5] = load_lcurves(R, repo_root, condition);
if ~isempty(T5)
    pr('Fig. 5: L-curve values from %s\n', src5);
    ttl = {'(a) common weight: \alpha_D = \alpha_S = \alpha', ...
           sprintf('(b) \\alpha_D sweep (\\alpha_S = 10^{%d})', round(log10(T5{2}(1, 3)))), ...
           sprintf('(c) \\alpha_S sweep (\\alpha_D = 10^{%d})', round(log10(T5{3}(1, 2))))};
    yl5 = {'H^1 seminorm of scaled D and S', 'H^1 seminorm of scaled D', 'H^1 seminorm of scaled S'};
    acol = [2 2 3];  % column of the swept weight: (a) alpha (= alphaD), (b) alphaD, (c) alphaS
    figH = figure('Color', 'w', 'Position', [60 60 1500 480]);
    for p = 1:3
        T = T5{p};  a = T(:, acol(p));  mis = T(:, 4);  sn = T(:, 5);
        subplot(1, 3, p);
        loglog(mis, sn, 'k-o', 'LineWidth', 1.3, 'MarkerSize', 7, 'MarkerFaceColor', 'w'); hold on;
        for k = 1:numel(a)
            text(mis(k)*1.03, sn(k)*1.25, sprintf('10^{%d}', round(log10(a(k)))), 'FontSize', 9);
        end
        kc = lc_corner(a, mis, sn);                      % maximum curvature (diamond)
        if ~isnan(kc), loglog(mis(kc), sn(kc), 'd', 'MarkerSize', 15, 'Color', [0.45 0.45 0.45], 'LineWidth', 1.3); end
        ks = find(abs(a - sel5(p)) <= 1e-6*sel5(p), 1);   % selected value (square)
        if ~isempty(ks), loglog(mis(ks), sn(ks), 's', 'MarkerSize', 17, 'Color', 'r', 'LineWidth', 1.6); end
        hold off; grid on; box on;
        xl = [min(mis)/1.15, max(mis)*1.35];  xt = [0.05 0.07 0.1 0.15 0.2 0.3 0.5 0.7];
        xt = xt(xt >= xl(1) & xt <= xl(2));
        xlim(xl);  ylim([min(sn)/2.5, max(sn)*2.5]);
        set(gca, 'XTick', xt, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), xt, 'UniformOutput', false));
        xlabel('data misfit'); ylabel(yl5{p}, 'Interpreter', 'tex');
        title(ttl{p}, 'Interpreter', 'tex', 'FontWeight', 'normal');
        set(gca, 'FontSize', fs - 1);
        if ~isnan(kc), pr('  panel %c: maximum curvature at %.0e, selected %.0e\n', 'a' + p - 1, a(kc), sel5(p)); end
    end
    save_figure_safe(figH, 'Fig05_L_curves.png', 300);  done{end+1} = 'Fig05_L_curves.png';
else
    skipped{end+1} = sprintf(['Fig. 5 (no L-curve values: run step 2 with alpha_mode = ''lcurve2'' and ', ...
                              'lcurve_only = true, or use condition T28 with data/T28/lcurve_article.csv)']);
end

%% Figs. 6-7: iteration histories of the four CG methods
mcol = {[1 0 0], [0 0 1], [0 0.75 0], [1 0 1]};   % FR, PR, PB, DY
mmrk = {'o', 's', 'd', '^'};
hist_f = {'J_history', 'grad_norm_history'; 'forward_error_history', 'inverse_error_history'};
hist_y = {'total cost J', '||\nablaJ||'; 'forward error', 'inverse error'};
hist_log = [true true; false true];   % semilogarithmic y axis as in the article
hist_k0  = [0 0; 0 1];                % first iteration shown (inverse error from k = 1)
hist_png = {'Fig06_convergence.png', 'Fig07_errors.png'};
for f = 1:2
    figH = figure('Color', 'w', 'Position', [60 60 1400 520]);
    for p = 1:2
        subplot(1, 2, p); hold on;
        h = gobjects(1, numel(cg_methods));  kmax = 0;
        for i = 1:numel(cg_methods)
            v = R.all_results{i}.(hist_f{f, p})(:).';  k = 0:numel(v) - 1;
            keep = k >= hist_k0(f, p) & isfinite(v);
            j = find(strcmpi({'FR', 'PR', 'PB', 'DY'}, cg_methods{i}), 1);  if isempty(j), j = i; end
            h(i) = plot_lm(k(keep), v(keep), mcol{j}, '-', mmrk{j}, 15);
            kmax = max(kmax, k(end));
        end
        hold off; grid on; box on; xlim([0 kmax]);
        if hist_log(f, p), set(gca, 'YScale', 'log'); end
        xlabel('iteration'); ylabel(hist_y{f, p}, 'Interpreter', 'tex');
        title(char('(a)' + (p - 1)*[0 1 0]), 'FontWeight', 'normal');
        legend(h, cg_methods, 'Location', 'northeast');
        set(gca, 'FontSize', fs + 1);
    end
    save_figure_safe(figH, hist_png{f}, 300);  done{end+1} = hist_png{f}; %#ok<SAGROW>
end

%% Figs. 8(a), 10(a): recovered D and S, 3D surfaces on a common scale
lab = {'D', 'S'};  F = {D, S};  unitF = {'mm^2/s', 'mol/(L s)'};
png3 = {sprintf('Fig08a_D_3D_%s.png', fig_method), sprintf('Fig10a_S_3D_%s.png', fig_method)};
for q = 1:2
    allv = [F{q}{1}(:); F{q}{2}(:); F{q}{3}(:)];  lims = [min(allv) max(allv)];
    if lims(2) <= lims(1), lims = lims + [-1 1]*eps; end
    figH = figure('Color', 'w', 'Position', [60 60 1500 500]);
    for s = 1:3
        ax = axes('Position', [0.04 + (s - 1)*0.30, 0.12, 0.25, 0.78]);
        surf(X, Y, F{q}{s}, 'EdgeColor', 'none', 'FaceColor', 'interp');
        colormap(ax, parula);  caxis(lims);  zlim(lims);
        xlim([min(x) max(x)]);  ylim([min(y) max(y)]);
        view(45, 30);  lighting gouraud;  camlight('headlight');  grid on;
        xlabel('x (mm)');  ylabel('y (mm)');
        zlabel(sprintf('%s_%d (%s)', lab{q}, s, unitF{q}), 'Interpreter', 'tex');
        title(sprintf('%s_%d (%s)', lab{q}, s, fig_method), 'Interpreter', 'tex');
        set(ax, 'FontSize', fs - 1);
    end
    cb = colorbar(ax);  set(cb, 'Position', [0.95 0.15 0.012 0.7]);
    cbar_label(cb, sprintf('%s (%s)', lab{q}, unitF{q}));
    save_figure_safe(figH, png3{q}, 300);  done{end+1} = png3{q}; %#ok<SAGROW>
end

%% Figs. 8(b), 10(b): recovered D and S, transverse profiles at y = 1, 4, 8 mm
pngl = {sprintf('Fig08b_D_lines_%s.png', fig_method), sprintf('Fig10b_S_lines_%s.png', fig_method)};
ylabl = {'D (10^{-3} mm^2/s)', 'S (10^{-2} mol/(L s))'};  scl = [scD scS];
for q = 1:2
    figH = figure('Color', 'w', 'Position', [60 60 1500 450]);
    for s = 1:3
        subplot(1, 3, s); hold on;
        h = gobjects(1, 3);
        for k = 1:3
            h(k) = plot(x, F{q}{s}(iy(k), :)*scl(q), lstyle{k}, 'Color', cols{s}, 'LineWidth', 1.8);
        end
        leg = ylab;
        if q == 1
            h(4) = plot([-0.5 0.5], D_lit(s)*scD*[1 1], 'k-.', 'LineWidth', 1.1);  leg{4} = 'reference';
        else
            plot([-0.5 0.5], [0 0], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
        end
        hold off; grid on; box on; xlim([-0.5 0.5]); set(gca, 'XTick', -0.5:0.25:0.5, 'FontSize', fs);
        xlabel('x (mm)'); title(species{s}, 'Interpreter', 'tex', 'FontWeight', 'normal');
        if s == 1, ylabel(ylabl{q}, 'Interpreter', 'tex'); end
        legend(h, leg, 'Location', 'best');
    end
    save_figure_safe(figH, pngl{q}, 300);  done{end+1} = pngl{q}; %#ok<SAGROW>
end

%% Fig. 9: simultaneous (fig_method) vs sequential recovery (D fixed at the reference values)
% line style: method (solid Sim., dashed Seq.); marker: y station (readable in black and white).
% Two short legends inside the axes: method in the first panel, y station in the last panel.
if isfield(R, 'seq_result') && ~isempty(R.seq_result)
    sq = R.seq_result;
    figH = figure('Color', 'w', 'Position', [60 60 1500 460]);
    for s = 1:3
        ax = axes('Position', [0.06 + (s - 1)*0.325, 0.14, 0.27, 0.76]); hold on;
        Sq = seq_field(sq, s);
        for k = 1:3
            plot_lm(x, S{s}(iy(k), :)*scS, ycols{k}, '-',  ymark{k}, 8);
            plot_lm(x, Sq(iy(k), :)*scS,   ycols{k}, '--', ymark{k}, 8);
        end
        plot([-0.5 0.5], [0 0], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
        if s == 1        % legend for the method (line style)
            hm = [plot(NaN, NaN, '-', 'Color', 'k', 'LineWidth', 1.6), plot(NaN, NaN, '--', 'Color', 'k', 'LineWidth', 1.6)];
            legend(ax, hm, {'Sim.', 'Seq.'}, 'Location', 'best', 'FontSize', fs - 1);
        elseif s == 3    % legend for the y station (marker)
            hy = gobjects(1, 3);
            for k = 1:3
                hy(k) = plot(NaN, NaN, ymark{k}, 'Color', ycols{k}, 'MarkerSize', 6, 'LineWidth', 1.1);
            end
            legend(ax, hy, ylab, 'Location', 'best', 'FontSize', fs - 1);
        end
        hold off; grid on; box on; xlim([-0.5 0.5]); set(gca, 'XTick', -0.5:0.25:0.5, 'FontSize', fs);
        xlabel('x (mm)'); title(species{s}, 'Interpreter', 'tex', 'FontWeight', 'normal');
        if s == 1, ylabel('S (10^{-2} mol/(L s))', 'Interpreter', 'tex'); end
    end
    save_figure_safe(figH, sprintf('Fig09_seq_S_%s.png', fig_method), 300);
    done{end+1} = sprintf('Fig09_seq_S_%s.png', fig_method);
    pr('Sequential recovery: data misfit %s %.4e, sequential %.4e (ratio %.2f)\n', fig_method, ...
       res.J_data_history(end), sq.J_data_history(end), sq.J_data_history(end)/res.J_data_history(end));
else
    skipped{end+1} = 'Fig. 9 (no seq_result: run step 2 with run_sequential_control = true)';
end

%% Numbers for Table 2 (D) and Section 6.4 (S)
bulk_mask = {X < -x_bulk, X > x_bulk, abs(X) > x_bulk};   % outside the reaction zone (Table 2)
pr('--- Table 2: recovered D (mm^2/s), method %s ---\n', fig_method);
pr('bulk: side where the species is present, x < -%.2f (H2SO4), x > %.2f (NaOH), |x| > %.2f mm (Na2SO4);\n', ...
   x_bulk, x_bulk, x_bulk);
pr('interface: median |x| < %.2f mm, max |x| < %.2f mm\n', x_int, x_max_int);
pr('%-7s | %-24s | %-24s | %-10s | %s\n', 'species', 'bulk median [p10, p90]', 'interface median [max]', 'int/bulk', 'reference');
for s = 1:3
    bulk = bulk_mask{s};
    intf = abs(X) < x_int;
    db = D{s}(bulk);  di = D{s}(intf);  dmax = max(D{s}(abs(X) < x_max_int));
    if isempty(db)
        pr('%-7s | %-24s | %.2e [%.2e]      | %-10s | %.2e\n', names{s}, 'no bulk points', ...
           median(di), dmax, '-', D_lit(s));
    else
        pr('%-7s | %.2e [%.2e, %.2e] | %.2e [%.2e]      | %6.2f     | %.2e\n', names{s}, median(db), ...
           prctile_(db, 10), prctile_(db, 90), median(di), dmax, median(di)/median(db), D_lit(s));
    end
end

pr('\n--- Section 6.4: source terms, method %s ---\n', fig_method);
I = zeros(1, 3);
for s = 1:3
    I(s) = trapz(y, trapz(x, S{s}, 2));
    Phi = trapz(x, S{s}, 2);
    pr('%-7s S range [%.3e, %.3e] mol/(L s); integral over domain %.4e mol/(L s) mm^2; ', names{s}, ...
       min(S{s}(:)), max(S{s}(:)), I(s));
    pr('Phi(y) at y = 1/4/8 mm: %.3e / %.3e / %.3e mol/(L s) mm\n', ...
       interp1(y, Phi, 1), interp1(y, Phi, 4), interp1(y, Phi, 8));
end
pr('Integrated ratio S1 : S2 : S3 = 1 : %.2f : %.2f   (stoichiometry 1 : 2 : -1)\n', I(2)/I(1), I(3)/I(1));
for k = 1:3
    Phi = cellfun(@(f) interp1(y, trapz(x, f, 2), y_lines(k)), S);
    pr('  at y = %d mm (width-integrated): 1 : %.2f : %.2f\n', y_lines(k), Phi(2)/Phi(1), Phi(3)/Phi(1));
end
neg2 = mean(S{2}(:) > 0);
pr('Fraction of the domain where S(NaOH) > 0: %.1f%%\n', 100*neg2);
pr('Data misfit %.4e, final J %.4e\n', res.J_data_history(end), res.final_cost);

%% Numbers for Section 6.2: comparison of the four CG methods
% Same definitions as Table 2 (bulk: outside the reaction zone) and Section 6.4
% (integrals over the domain); interfacial region |x| < x_max_int. The local difference
% of a field from fig_method is the relative L2 norm ||F_m - F_ref|| / ||F_ref|| over |x| < x_max_int.
pr('\n--- Section 6.2: comparison of the CG methods (reference %s) ---\n', fig_method);
nm   = numel(cg_methods);
intf = abs(X) < x_max_int;
I3   = zeros(1, nm);  rat = zeros(nm, 2);  D3b = nan(1, nm);
pr('%-4s | %-8s | %-8s | %-22s | %-14s | %-16s\n', 'meth', 'J', 'misfit', 'ratio S1:S2:S3 (int.)', ...
   'int. S3', 'D3 bulk median');
for i = 1:nm
    ri = R.all_results{i};
    Si = {ri.S1_phys, ri.S2_phys, ri.S3_phys};
    Ii = cellfun(@(f) trapz(y, trapz(x, f, 2)), Si);
    I3(i) = Ii(3);  rat(i, :) = [Ii(2) Ii(3)]/Ii(1);
    D3b(i) = median(ri.D3_phys(bulk_mask{3}));
    pr('%-4s | %8.4f | %8.4f | 1 : %5.2f : %5.2f      | %.4e    | %.3e mm^2/s\n', cg_methods{i}, ...
       ri.final_cost, ri.J_data_history(end), rat(i, 1), rat(i, 2), I3(i), D3b(i));
end
pr('Integrated S3: %.4e to %.4e; max deviation from the mean %.1f%%\n', min(I3), max(I3), ...
   100*max(abs(I3 - mean(I3)))/mean(I3));
pr('D3 bulk median: %.2f to %.2f x 10^-3 mm^2/s\n', 1e3*min(D3b), 1e3*max(D3b));
pr('Ratio S2/S1: %.2f to %.2f; S3/S1: %.2f to %.2f\n', min(rat(:, 1)), max(rat(:, 1)), ...
   min(rat(:, 2)), max(rat(:, 2)));
pr('Relative L2 difference from %s in |x| < %.2f mm:\n', fig_method, x_max_int);
pr('%-4s | %6s %6s %6s | %6s %6s %6s\n', 'meth', 'D1', 'D2', 'D3', 'S1', 'S2', 'S3');
F0 = {res.D1_phys, res.D2_phys, res.D3_phys, res.S1_phys, res.S2_phys, res.S3_phys};
dall = [];
for i = 1:nm
    if i == m, continue; end
    ri = R.all_results{i};
    Fi = {ri.D1_phys, ri.D2_phys, ri.D3_phys, ri.S1_phys, ri.S2_phys, ri.S3_phys};
    dl = cellfun(@(a, b) norm(b(intf) - a(intf))/norm(a(intf)), F0, Fi);
    dall = [dall dl]; %#ok<AGROW>
    pr('%-4s | %6.2f %6.2f %6.2f | %6.2f %6.2f %6.2f\n', cg_methods{i}, dl);
end
if ~isempty(dall)
    pr('Range of the local differences: %.0f%% to %.0f%%\n', 100*min(dall), 100*max(dall));
end

%% Fig. 11: Monte Carlo standard deviations (if available)
if isfield(R, 'Ds') && ~isempty(R.Ds)
    figH = figure('Color', 'w', 'Position', [60 60 1500 800]);
    unit11 = {'mm^2/s', 'mol/(L s)'};  MC = {R.Ds, R.Ss};
    for q = 1:2
        for s = 1:3
            sd = std(MC{q}(:, :, s, :), 0, 4);
            subplot(2, 3, 3*(q - 1) + s); imagesc(x, y, sd); axis xy; colorbar;
            ylim([0 max(y)]); set(gca, 'FontSize', fs);
            xlabel('x (mm)'); if s == 1, ylabel('y (mm)'); end
            title(sprintf('std of %s_%d, %s (%s)', lab{q}, s, species{s}, unit11{q}), 'Interpreter', 'tex');
        end
    end
    save_figure_safe(figH, 'Fig11_MC_std.png', 300);  done{end+1} = 'Fig11_MC_std.png';
else
    skipped{end+1} = 'Fig. 11 (no Monte Carlo fields: run step 2 with run_items.mc = true)';
end
fclose(fid);
fprintf('\nNumbers written to step3_article_numbers.txt\n');

%% Figs. 12-13: D and S for the temperature series (20, 28, 37 degC) at y = t_y
cd(start_dir);
have = false(1, numel(t_conds));
for c = 1:numel(t_conds)
    have(c) = exist(fullfile(repo_root, 'results', t_conds{c}, 'step2_cg_results.mat'), 'file') == 2 && ...
              exist(fullfile(repo_root, 'results', t_conds{c}, 'step1_initial_estimates.mat'), 'file') == 2;
end
if any(have)
    TS = struct('lab', {}, 'col', {}, 'ls', {}, 'mk', {}, 'x', {}, 'y', {}, 'D', {}, 'S', {}, 'Dlit', {});
    for c = find(have)
        Rk = load(fullfile(repo_root, 'results', t_conds{c}, 'step2_cg_results.mat'), 'all_results', 'cg_methods', 'x', 'y');
        Lk = load(fullfile(repo_root, 'results', t_conds{c}, 'step1_initial_estimates.mat'), 'D_lit_array', 'T_C');
        ik = find(strcmpi(cellstr(Rk.cg_methods), fig_method), 1);
        rk = Rk.all_results{ik};
        TS(end+1) = struct('lab', sprintf('%d \\circC', round(double(Lk.T_C))), 'col', t_cols{c}, ...
                           'ls', t_ls{c}, 'mk', t_mark{c}, 'x', Rk.x(:).', 'y', Rk.y(:), ...
                           'D', {{rk.D1_phys, rk.D2_phys, rk.D3_phys}}, ...
                           'S', {{rk.S1_phys, rk.S2_phys, rk.S3_phys}}, 'Dlit', Lk.D_lit_array); %#ok<SAGROW>
    end
    q12 = {'D', 'S'};  sc12 = [scD scS];
    png12 = {'Fig12_D_temperature.png', 'Fig13_S_temperature.png'};
    for f = 1:2
        figH = figure('Color', 'w', 'Position', [60 60 1500 480]);
        for s = 1:3
            axes('Position', [0.07 + (s - 1)*0.315, 0.15, 0.26, 0.72]); hold on;
            h = gobjects(1, numel(TS));
            for c = 1:numel(TS)
                jy = nearest_index(TS(c).y, t_y);
                h(c) = plot_lm(TS(c).x, TS(c).(q12{f}){s}(jy, :)*sc12(f), TS(c).col, TS(c).ls, TS(c).mk, 8);
                if f == 1
                    plot([-0.5 0.5], TS(c).Dlit(s)*scD*[1 1], ':', 'Color', TS(c).col, 'LineWidth', 1.1);
                end
            end
            if f == 2, plot([-0.5 0.5], [0 0], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8); end
            hold off; grid on; box on; xlim([-0.5 0.5]); set(gca, 'XTick', -0.5:0.25:0.5, 'FontSize', fs);
            title(species{s}, 'Interpreter', 'tex', 'FontWeight', 'normal');
            xlabel('x (mm)');
            if s == 1
                ylabel(sprintf('%s (10^{%d} %s)', q12{f}, -round(log10(sc12(f))), unitF{f}), 'Interpreter', 'tex');
                legend(h, {TS.lab}, 'Location', 'northwest');
            end
        end
        save_figure_safe(figH, fullfile(repo_root, 'results', png12{f}), 300);  done{end+1} = png12{f}; %#ok<SAGROW>
    end
    if ~all(have)
        fprintf(['Figs. 12-13 drawn with %s only. Run steps 1 and 2 with condition = %s ', ...
                 '(default settings of step 2) to add the missing temperatures.\n'], ...
                strjoin(t_conds(have), ', '), strjoin(t_conds(~have), ' and '));
    end
else
    skipped{end+1} = 'Figs. 12-13 (run steps 1 and 2 for T20, T28 and T37 first)';
end
cd(start_dir);

fprintf('\n=== Step 3: %d figure files written ===\n', numel(done));
fprintf('  %s\n', done{:});
if ~isempty(skipped)
    fprintf('Not drawn:\n');  fprintf('  %s\n', skipped{:});
end

%% ---------------- helpers ----------------
function i = nearest_index(v, v0)
    [~, i] = min(abs(v(:) - v0));
end

function B = medfilt3x3(A)
% 3 x 3 median filter with symmetric padding (as medfilt2(A, [3 3], 'symmetric'), base MATLAB)
    P = A([1 1:end end], [1 1:end end]);
    [n1, n2] = size(A);  Z = zeros(n1, n2, 9);  k = 0;
    for i = 0:2
        for j = 0:2
            k = k + 1;  Z(:, :, k) = P(1 + i:n1 + i, 1 + j:n2 + j);
        end
    end
    B = median(Z, 3);
end

function cbar_label(cb, str)
% colour-bar label (MATLAB colorbar object; axes-based colour bars as fallback)
    try, cb.Label.String = str; cb.Label.Interpreter = 'tex';
    catch, try, xlabel(cb, str); catch, end
    end
end

function cbar_ticks(cb, t)
    try, set(cb, 'Ticks', t); catch, try, set(cb, 'XTick', t); catch, end
    end
end

function h = plot_lm(x, f, col, ls, mk, every)
% line with a marker every EVERY points; returns a handle that shows both in the legend
    plot(x, f, ls, 'Color', col, 'LineWidth', 1.6);
    idx = 1:every:numel(x);
    plot(x(idx), f(idx), mk, 'Color', col, 'MarkerSize', 6, 'LineWidth', 1.1);
    h = plot(NaN, NaN, [ls mk], 'Color', col, 'LineWidth', 1.6, 'MarkerSize', 6);
end

function Sq = seq_field(sq, s)
% S of the sequential recovery (field names differ between versions of step 2)
    if isfield(sq, sprintf('S%d_phys', s)), Sq = sq.(sprintf('S%d_phys', s));
    else, Sq = sq.(sprintf('S%d', s)); end
end

function [T, sel, src] = load_lcurves(R, repo_root, condition)
% L-curve tables [panel, alphaD, alphaS, misfit, seminorm] for Fig. 5 and the selected weights
    T = {};  sel = [NaN NaN NaN];  src = '';
    LC = [];  LC2 = [];
    f1 = fullfile(repo_root, 'results', condition, 'step2_lcurve_results.mat');
    if exist(f1, 'file') == 2
        M = load(f1);  LC = M.LC;  LC2 = M.LC2;  src = f1;
    elseif isfield(R, 'LC') && ~isempty(R.LC) && isfield(R, 'LC2') && ~isempty(R.LC2)
        LC = R.LC;  LC2 = R.LC2;  src = 'LC/LC2 in step2_cg_results.mat';
    end
    if ~isempty(LC) && ~isempty(LC2)
        a = LC.table(:, 1);
        T{1} = [ones(size(a)), a, a, LC.table(:, 2), LC.table(:, 3)];
        sw = LC2.sweeps;  if ~iscell(sw), sw = num2cell(sw); end
        iD = find(cellfun(@(w) strcmpi(w.which, 'D'), sw), 1);
        iS = find(cellfun(@(w) strcmpi(w.which, 'S'), sw), 1);
        wD = sw{iD};  wS = sw{iS};
        T{2} = [2*ones(size(wD.table, 1), 1), wD.table(:, 1), wD.fixed*ones(size(wD.table, 1), 1), wD.table(:, 2:3)];
        T{3} = [3*ones(size(wS.table, 1), 1), wS.fixed*ones(size(wS.table, 1), 1), wS.table(:, 1), wS.table(:, 2:3)];
        sel = [LC.alpha_discrepancy, R.alphaD, R.alphaS];
        return
    end
    f3 = fullfile(repo_root, 'data', condition, 'lcurve_article.csv');
    if exist(f3, 'file') == 2
        A = dlmread(f3, ',', 1, 0);
        for p = 1:3, T{p} = A(A(:, 1) == p, 1:5); end
        flag = A(:, 6);
        sel(1) = A(A(:, 1) == 1 & flag == 1, 2);
        sel(2) = A(A(:, 1) == 2 & flag == 1, 2);
        sel(3) = A(A(:, 1) == 3 & flag == 1, 3);
        src = [f3, ' (values of the article run)'];
    end
end

function k = lc_corner(a, mis, sn)
% index of the maximum curvature of (log misfit, log seminorm) vs log alpha, interior points
% only (same rule as lcurve_choice in step 2)
    k = NaN;
    if numel(a) < 3, return; end
    la = log10(a(:));  lr = log10(mis(:));  ls = log10(max(sn(:), realmin));
    r1 = gradient(lr, la);  s1 = gradient(ls, la);
    r2 = gradient(r1, la);  s2 = gradient(s1, la);
    kap = (r1.*s2 - r2.*s1)./max((r1.^2 + s1.^2).^1.5, eps);
    [~, i] = max(kap(2:end-1));
    k = i + 1;
end

function p = prctile_(v, q)
% percentile without the Statistics Toolbox
    v = sort(v(:));  n = numel(v);
    p = v(max(1, min(n, round(q/100*(n - 1)) + 1)));
end

function save_figure_safe(figHandle, fname, dpi)
% exportgraphics if available, otherwise print; never throws
    drawnow;  ok = false;
    if exist('exportgraphics', 'file') == 2
        try, exportgraphics(figHandle, fname, 'Resolution', dpi); ok = true; catch, end
    end
    if ~ok
        try, print(figHandle, fname, '-dpng', sprintf('-r%d', dpi)); ok = true; catch, end
    end
    if ~ok
        try, print(figHandle, fname, '-dpng', '-r150'); ok = true; catch, end
    end
    if ok, fprintf('Saved %s\n', fname); else, warning('Could not save %s', fname); end
end
