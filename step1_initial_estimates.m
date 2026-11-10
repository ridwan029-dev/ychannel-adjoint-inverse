%STEP1_INITIAL_ESTIMATES  Profile fitting and direct initial estimates of D and S.
%
%   Step 1 of 3. Reads the measured NIR concentration maps of H2SO4, NaOH and
%   Na2SO4 of one data set, centres the profiles at the Na2SO4 maximum, fits
%   smooth erfc/Gaussian profiles by least squares, computes the spatial
%   derivatives, the initial source terms from the width-integrated balance
%   (Eqs. (17)-(19)) and the direct estimate of D from the steady
%   advection-diffusion-reaction equation (Eqs. (20)-(21)), and checks and
%   plots the results.
%
%   Input    data/<condition>/H2SO4.tif, NaOH.tif, Na2SO4.tif
%            measured concentration maps (mol/L); <condition> is set below
%   Output   in results/<condition>/:
%            step1_initial_estimates.mat  fitted fields, initial D and S (input of step 2)
%            step1_intermediate.mat       intermediate fields and derivatives
%            map_concentration_medfilt.mat  measured maps after a 3 x 3 median filter (Fig. 2(b))
%            *.png                        figures of the direct estimate
%   Results  Figs. 2 and 4, Section 6.1 (condition T28)
%   Usage    Set the condition below, set the MATLAB current folder to the
%            repository root and run
%              >> step1_initial_estimates
%   Requires MATLAB R2020a or newer, Image Processing Toolbox (imgaussfilt),
%            Curve Fitting Toolbox (smooth)
%   Runtime  A few minutes
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

clear all; close all; clc;

%% ========================================================================
%% SETTINGS: DATA SET
%% ========================================================================
% NIR concentration maps measured with the setup of Uema et al. (2021).
% T28 is the data set of Sections 6.1-6.5; T20, T28 and T37 are used in Section 6.6.
condition = 'T28';
%            tag     v (mm/s)  T (degC)
cond_table = {'T20',  0.15,     20;     % disk heater
              'T28',  0.15,     28;
              'T37',  0.15,     37};
ic = find(strcmpi(cond_table(:, 1), condition), 1);
assert(~isempty(ic), 'Unknown condition %s.', condition);
v_cond    = cond_table{ic, 2};          % mean velocity, mm/s
T_C       = cond_table{ic, 3};          % temperature, degC
pixel_mm  = 0.03;                       % pixel size in x and y (30 um pixel pitch, 1x lens)

start_dir = pwd;
repo_root = fileparts(mfilename('fullpath'));
if isempty(repo_root), repo_root = start_dir; end
image_dir = fullfile(repo_root, 'data', condition);
out_dir   = fullfile(repo_root, 'results', condition);   % all results of this data set
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
cd(out_dir);
fprintf('Results folder: %s\n', out_dir);

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  STEP 1: DATA LOADING, ERFC PROFILES & DIFFUSION ESTIMATION   ║\n');
fprintf('║  Profiles: Smooth ERFC | Diffusion: Direct PDE iterative      ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

%% ========================================================================
%% SECTION 1: REFERENCE VALUES OF D
%% ========================================================================
fprintf('=== LITERATURE REFERENCE VALUES (%d°C) ===\n', T_C);

% Values at 30 °C scaled to T_C with the Stokes-Einstein relation, D ~ T/eta(T)
eta_w = @(T) 2.414e-5 * 10.^(247.8 ./ (T + 273.15 - 140));   % water viscosity, Pa s (Vogel)
D30   = [1.85e-3, 2.13e-3, 0.85e-3];                           % mm^2/s at 30 °C
D_T   = D30 * ((T_C + 273.15)/303.15) * (eta_w(30)/eta_w(T_C));
D_lit = struct();
D_lit.H2SO4  = D_T(1);
D_lit.NaOH   = D_T(2);
D_lit.Na2SO4 = D_T(3);

fprintf('  D(H₂SO₄)  = %.2e mm²/s\n', D_lit.H2SO4);
fprintf('  D(NaOH)   = %.2e mm²/s\n', D_lit.NaOH);
fprintf('  D(Na₂SO₄) = %.2e mm²/s\n\n', D_lit.Na2SO4);

D_lit_array = [D_lit.H2SO4, D_lit.NaOH, D_lit.Na2SO4];

%% ========================================================================
%% SECTION 2: NON-DIMENSIONALIZATION WITH VERIFICATION
%% ========================================================================
fprintf('=== NON-DIMENSIONALIZATION PARAMETERS ===\n');

ndV  = 0.1;    % Velocity scale:      v* = 1 × 10⁻¹ (mm/s)
ndL  = 0.01;   % Length scale:        L* = 1 × 10⁻² (mm)
ndD  = 0.001;  % Diffusion scale:     D* = 1 × 10⁻³ (mm²/s)
ndC  = 1;      % Concentration scale: C* = 1 (mol/L)
ndS  = ndV*ndC/ndL;  % Source-term scale:   S* = v C / L = 10 mol/(L s)
use_full_laplacian = false;  % false: d2C/dx2 only, as in Eq. (2); true: d2C/dx2 + d2C/dy2
dpi_ = pixel_mm;   % pixel size (mm)

sig_nd        = dpi_ / ndL;          % σ* = 3
v_dimensional = v_cond;              % mean velocity (mm/s)
v_nd          = v_dimensional / ndV;

% Péclet number
Pe_typical = v_dimensional * 0.4 / D_lit.H2SO4;
fprintf('  Péclet number (typical): Pe = %.1f\n', Pe_typical);

% Consistency check
fprintf('\n  Non-dim consistency check:\n');
fprintf('    ndV × ndL = %.4f\n', ndV * ndL);
fprintf('    ndD       = %.4f\n', ndD);
if abs(ndV * ndL - ndD) < 1e-10
    fprintf('    ✓ Scaling is consistent\n\n');
else
    fprintf('    ⚠ Warning: Scaling inconsistency detected!\n');
    fprintf('    Consider: ndD = ndV × ndL = %.4f\n\n', ndV * ndL);
end

fprintf('  Scale factors:\n');
fprintf('    ndV = %.3f mm/s\n',      ndV);
fprintf('    ndL = %.3f mm\n',        ndL);
fprintf('    ndD = %.6f mm²/s\n',     ndD);
fprintf('    ndS = %.3f mol/L·s\n',   ndS);
fprintf('    ndC = %.3f mol/L\n\n',   ndC);

%% ========================================================================
%% SECTION 3: LOAD AND PREPARE DATA
%% ========================================================================
fprintf('=== LOADING EXPERIMENTAL DATA (%s: v = %.2f mm/s, T = %d °C) ===\n', condition, v_cond, T_C);
species_tags_file = {'H2SO4', 'NaOH', 'Na2SO4'};

if ~exist(fullfile(image_dir, 'H2SO4.tif'), 'file')
    cd(start_dir);
    error('Measured concentration maps not found in %s (expected H2SO4.tif, NaOH.tif, Na2SO4.tif).', image_dir);
end
% After fliplr, image rows run across the channel (x) and columns 51:end along it (y).
% F{s} is Ny x 256: one row per y position, one column per image row.
F = cell(1, 3);
for s = 1:3
    im = double(imread(fullfile(image_dir, [species_tags_file{s} '.tif']))) / ndC;
    im = fliplr(im);
    F{s} = im(:, 51:end)';
end
fprintf('Images loaded from: %s\n', image_dir);

%% ========================================================================
%% SECTION 4: DEFINE COMPUTATIONAL GRID
%% ========================================================================
fprintf('\n=== DEFINING COMPUTATIONAL GRID ===\n');

[ny, nx] = deal(size(F{1}, 1), 81);
x_min = -0.5;   % mm, relative to the interface
x_max =  0.5;   % mm
x_phys_mm = linspace(x_min, x_max, nx);
y_phys_mm = (0:ny-1) * pixel_mm;        % one grid row per image column
y_max = y_phys_mm(end);

% Centre-position adjustment (Carandang et al., 2024, Sec. 3.1): at every y the profiles are
% translated so that the Na2SO4 maximum lies at x = 0; the same shift is applied to H2SO4 and NaOH.
r_img    = 1:size(F{1}, 2);             % image row index across the channel
r_search = 91:180;                       % image rows that contain the interface
mu_row   = zeros(ny, 1);
for j = 1:ny
    jj = max(1, j-2):min(ny, j+2);       % average over 5 neighbouring y positions
    p  = mean(F{3}(jj, r_search), 1);
    [~, k] = max(p);
    mu_row(j) = r_search(k);
    % Gaussian fit a*exp(-(r - m)^2/(2 w^2)) + b around the maximum (no toolbox needed)
    kk = max(1, k-20):min(numel(p), k+20);
    rr = r_search(kk);  pp = p(kk);
    sse = @(q) sum((q(1)*exp(-(rr - q(2)).^2/(2*q(3)^2)) + q(4) - pp).^2);
    q = fminsearch(sse, [max(pp) - min(pp), r_search(k), 5, min(pp)], optimset('Display', 'off'));
    if abs(q(2) - r_search(k)) < 10, mu_row(j) = q(2); end
end
mu_row = movmean(mu_row, 15, 'Endpoints', 'shrink');   % smooth shift along y
fprintf('Interface (image row): %.1f at y = 0, %.1f at y = 4 mm, %.1f at y = %.1f mm\n', ...
        mu_row(1), mu_row(round(4/pixel_mm) + 1), mu_row(end), y_max);

C1_obs = zeros(ny, nx);
C2_obs = zeros(ny, nx);
C3_obs = zeros(ny, nx);
for j = 1:ny
    r = mu_row(j) + x_phys_mm / pixel_mm;
    C1_obs(j,:) = interp1(r_img, F{1}(j,:), r, 'linear');
    C2_obs(j,:) = interp1(r_img, F{2}(j,:), r, 'linear');
    C3_obs(j,:) = interp1(r_img, F{3}(j,:), r, 'linear');
end
x_map_mode = 'pixel, centred at each y';

% Thresholds of the validity mask of the direct estimate (Section 5)
lap_thr    = 1e-3 * (5/12)^2;            % |d2C*/dx*2| > 1.7e-4
grad_thr   = 1e-3 * (5/12);              % |dC*/dx*|   > 4.2e-4
D_upper_nd = 50;                         % upper bound of D* (0.05 mm^2/s)

% Positive constraint
C1_obs(C1_obs < 0) = 0;
C2_obs(C2_obs < 0) = 0;
C3_obs(C3_obs < 0) = 0;

% Mesh grids
x_nd = x_phys_mm / ndL;
y_nd = y_phys_mm / ndL;
[X_nd, Y_nd]     = meshgrid(x_nd, y_nd);
[X_phys, Y_phys] = meshgrid(x_phys_mm, y_phys_mm);

dx_nd   = x_nd(2) - x_nd(1);
dy_nd   = y_nd(2) - y_nd(1);
dx_phys = x_phys_mm(2) - x_phys_mm(1);
dy_phys = y_phys_mm(2) - y_phys_mm(1);

v_advection_nd = v_dimensional / ndV;

C_obs_all     = {C1_obs, C2_obs, C3_obs};
species_names = {'H₂SO₄', 'NaOH', 'Na₂SO₄'};
species_tags  = {'H2SO4',  'NaOH', 'Na2SO4'};

fprintf('Grid: %d × %d points\n',             ny, nx);
fprintf('Physical: x ∈ [%.2f, %.2f] mm (0 = interface), y ∈ [0, %.2f] mm\n', x_min, x_max, y_max);
fprintf('Resolution: dx = %.4f mm, dy = %.4f mm\n', dx_phys, dy_phys);
fprintf('Non-dimensional: v* = %.2f\n\n', v_advection_nd);

%% ========================================================================
%% SECTION 5: SMOOTH ERFC-BASED CONCENTRATION PROFILES
%% ========================================================================
fprintf('═══════════════════════════════════════════════════════════\n');
fprintf('=== SMOOTH ERFC CONCENTRATION PROFILES ===\n');
fprintf('═══════════════════════════════════════════════════════════\n\n');

profile_source = 'fit';   % least-squares fit to the measured maps
fit_info = struct();
fprintf('Least-squares fit of erfc/Gaussian profiles to the measured maps\n');
fprintf('(row-by-row fit, then each parameter linear in y)\n\n');
[C1_nd, P1, st1] = fit_profiles(C1_obs, x_nd, y_nd, 'erfc_left');
[C2_nd, P2, st2] = fit_profiles(C2_obs, x_nd, y_nd, 'erfc_right');
[C3_nd, P3, st3] = fit_profiles(C3_obs, x_nd, y_nd, 'gauss');
alpha1_interp = repmat(P1.a, 1, nx);  mu1_interp = repmat(P1.mu, 1, nx);  sigma1_interp = repmat(P1.s, 1, nx);
alpha2_interp = repmat(P2.a, 1, nx);  mu2_interp = repmat(P2.mu, 1, nx);  sigma2_interp = repmat(P2.s, 1, nx);
A3_interp     = repmat(P3.a, 1, nx);  mu3_interp = repmat(P3.mu, 1, nx);  sigma3_interp = repmat(P3.s, 1, nx);
x_center_S_nd = P3.mu;               % centre of the Gaussian S = fitted Na2SO4 peak position
names_fit = {'H2SO4 (erfc, left)', 'NaOH (erfc, right)', 'Na2SO4 (Gaussian)'};  Pall = {P1, P2, P3};  Sall = {st1, st2, st3};
for s = 1:3
    fprintf('  %-20s RMSE = %.4f mol/L, R^2 = %.3f (median row R^2 = %.3f)\n', names_fit{s}, Sall{s}.RMSE, Sall{s}.R2, Sall{s}.median_row_R2);
    fprintf('  %-20s a: %.3f -> %.3f mol/L,  mu: %.4f -> %.4f mm,  sigma: %.4f -> %.4f mm\n', '', ...
        Pall{s}.a(1), Pall{s}.a(end), Pall{s}.mu(1)*ndL, Pall{s}.mu(end)*ndL, Pall{s}.s(1)*ndL, Pall{s}.s(end)*ndL);
end
fprintf('  (NIR concentration resolution: 0.087 mol/L)\n\n');
fit_info.P = Pall;  fit_info.stats = Sall;
% quality of the fitted profiles (RMSE and R^2 against the measured maps)
C_tmp = {C1_nd, C2_nd, C3_nd};  Cm_tmp = {C1_obs, C2_obs, C3_obs};
for s = 1:3
    r_ = C_tmp{s}(:) - Cm_tmp{s}(:);
    fit_info.RMSE(s) = sqrt(mean(r_.^2));
    fit_info.R2(s) = 1 - sum(r_.^2)/sum((Cm_tmp{s}(:) - mean(Cm_tmp{s}(:))).^2);
end
fprintf('Profiles vs measured maps: RMSE = %.4f / %.4f / %.4f mol/L, R^2 = %.3f / %.3f / %.3f\n\n', fit_info.RMSE, fit_info.R2);
clear C_tmp Cm_tmp r_
C_nd_all = {C1_nd, C2_nd, C3_nd};

%% ========================================================================
%% SECTION 5b: PARAMETER SMOOTHNESS VERIFICATION
%% ========================================================================
fprintf('═══════════════════════════════════════════════════════════\n');
fprintf('=== PARAMETER SMOOTHNESS VERIFICATION ===\n');
fprintf('═══════════════════════════════════════════════════════════\n\n');

[~, x_center_idx] = min(abs(x_nd));

mu1_profile    = mu1_interp(:, x_center_idx);
sigma1_profile = sigma1_interp(:, x_center_idx);
mu2_profile    = mu2_interp(:, x_center_idx);
sigma2_profile = sigma2_interp(:, x_center_idx);
A3_profile     = A3_interp(:, x_center_idx);
sigma3_profile = sigma3_interp(:, x_center_idx);

dmu1_dy    = gradient(mu1_profile, dy_nd);
dsigma1_dy = gradient(sigma1_profile, dy_nd);
dmu2_dy    = gradient(mu2_profile, dy_nd);
dsigma2_dy = gradient(sigma2_profile, dy_nd);
dA3_dy     = gradient(A3_profile, dy_nd);
dsigma3_dy = gradient(sigma3_profile, dy_nd);

fprintf('Parameter gradient statistics (should be smooth):\n');
fprintf('  dμ₁/dy:  mean = %.4f, std = %.4f, max = %.4f\n', ...
        mean(dmu1_dy), std(dmu1_dy), max(abs(dmu1_dy)));
fprintf('  dσ₁/dy:  mean = %.4f, std = %.4f, max = %.4f\n', ...
        mean(dsigma1_dy), std(dsigma1_dy), max(abs(dsigma1_dy)));
fprintf('  dμ₂/dy:  mean = %.4f, std = %.4f, max = %.4f\n', ...
        mean(dmu2_dy), std(dmu2_dy), max(abs(dmu2_dy)));
fprintf('  dσ₂/dy:  mean = %.4f, std = %.4f, max = %.4f\n', ...
        mean(dsigma2_dy), std(dsigma2_dy), max(abs(dsigma2_dy)));
fprintf('  dA₃/dy:  mean = %.6f, std = %.6f, max = %.6f\n', ...
        mean(dA3_dy), std(dA3_dy), max(abs(dA3_dy)));
fprintf('  dσ₃/dy:  mean = %.6f, std = %.6f, max = %.6f\n', ...
        mean(dsigma3_dy), std(dsigma3_dy), max(abs(dsigma3_dy)));

CV_mu1    = std(dmu1_dy) / abs(mean(dmu1_dy));
CV_sigma1 = std(dsigma1_dy) / abs(mean(dsigma1_dy));

fprintf('\nSmoothness indicators (lower is smoother):\n');
fprintf('  CV(dμ₁/dy) = %.4f (should be < 0.1 for smooth)\n', CV_mu1);
fprintf('  CV(dσ₁/dy) = %.4f (should be < 0.1 for smooth)\n', CV_sigma1);

if CV_mu1 < 0.1 && CV_sigma1 < 0.1
    fprintf('  ✓ Parameters are SMOOTH\n\n');
else
    fprintf('  ⚠ Parameters may have discontinuities\n\n');
end

%% Na₂SO₄ concentration check
fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║          CRITICAL: Na₂SO₄ CONCENTRATION CHECK             ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

fprintf('  Min(C₃): %.6e mol/L\n', min(C3_nd(:)));
fprintf('  Max(C₃): %.6e mol/L\n', max(C3_nd(:)));
fprintf('  Mean(C₃): %.6e mol/L\n', mean(C3_nd(:)));

if max(C3_nd(:)) < 0.01
    fprintf('⚠⚠⚠ CRITICAL ERROR: C₃ is too small! ⚠⚠⚠\n\n');
elseif max(C3_nd(:)) > 1.0
    fprintf('⚠ WARNING: C₃ is too large (>1.0 mol/L)\n\n');
else
    fprintf('✓ C₃ concentration scale is reasonable\n\n');
end

%% ========================================================================
%% SECTION 6: SPATIAL DERIVATIVES (5-point stencil)
%% ========================================================================
fprintf('=== COMPUTING SPATIAL DERIVATIVES (5-point stencil) ===\n');

dCdx_nd_all   = cell(1,3);
dCdy_nd_all   = cell(1,3);
d2Cdxx_nd_all = cell(1,3);
d2Cdyy_nd_all = cell(1,3);

for s = 1:3
    fprintf('Computing non-dimensional derivatives for %s...\n', species_names{s});
    C_nd = C_nd_all{s};

    % Apply smoothing to reduce noise before differentiation
    C_nd_smooth = imgaussfilt(C_nd, 0.5);

    dCdx_nd  = zeros(ny, nx);
    dCdy_nd  = zeros(ny, nx);
    d2Cdxx_nd = zeros(ny, nx);
    d2Cdyy_nd = zeros(ny, nx);

    % ∂C*/∂x*: 5-point stencil interior, 3-point boundary
    if nx >= 5
        for i = 3:nx-2
            dCdx_nd(:,i) = (-C_nd_smooth(:,i+2) + 8*C_nd_smooth(:,i+1) ...
                            - 8*C_nd_smooth(:,i-1) + C_nd_smooth(:,i-2)) / (12*dx_nd);
        end
        dCdx_nd(:,1)    = (-3*C_nd_smooth(:,1) + 4*C_nd_smooth(:,2) - C_nd_smooth(:,3)) / (2*dx_nd);
        dCdx_nd(:,2)    = (-C_nd_smooth(:,1) + C_nd_smooth(:,3)) / (2*dx_nd);
        dCdx_nd(:,nx-1) = (C_nd_smooth(:,nx-2) - C_nd_smooth(:,nx)) / (2*dx_nd);
        dCdx_nd(:,nx)   = (C_nd_smooth(:,nx-2) - 4*C_nd_smooth(:,nx-1) + 3*C_nd_smooth(:,nx)) / (2*dx_nd);
    else
        for i = 2:nx-1
            dCdx_nd(:,i) = (C_nd_smooth(:,i+1) - C_nd_smooth(:,i-1)) / (2*dx_nd);
        end
        if nx > 1
            dCdx_nd(:,1)  = (C_nd_smooth(:,2)  - C_nd_smooth(:,1))    / dx_nd;
            dCdx_nd(:,nx) = (C_nd_smooth(:,nx) - C_nd_smooth(:,nx-1)) / dx_nd;
        end
    end

    % ∂C*/∂y*: 5-point stencil interior, 3-point boundary
    for j = 3:ny-2
        dCdy_nd(j,:) = (-C_nd_smooth(j+2,:) + 8*C_nd_smooth(j+1,:) ...
                        - 8*C_nd_smooth(j-1,:) + C_nd_smooth(j-2,:)) / (12*dy_nd);
    end
    dCdy_nd(1,:)    = (-3*C_nd_smooth(1,:) + 4*C_nd_smooth(2,:) - C_nd_smooth(3,:)) / (2*dy_nd);
    dCdy_nd(2,:)    = (-C_nd_smooth(1,:) + C_nd_smooth(3,:)) / (2*dy_nd);
    dCdy_nd(ny-1,:) = (C_nd_smooth(ny-2,:) - C_nd_smooth(ny,:)) / (2*dy_nd);
    dCdy_nd(ny,:)   = (C_nd_smooth(ny-2,:) - 4*C_nd_smooth(ny-1,:) + 3*C_nd_smooth(ny,:)) / (2*dy_nd);

    % ∂²C*/∂x*²: second derivatives (central difference)
    if nx >= 3
        for i = 2:nx-1
            d2Cdxx_nd(:,i) = (C_nd_smooth(:,i+1) - 2*C_nd_smooth(:,i) + C_nd_smooth(:,i-1)) / (dx_nd^2);
        end
        if nx >= 4
            d2Cdxx_nd(:,1)  = (2*C_nd_smooth(:,1)  - 5*C_nd_smooth(:,2)    + 4*C_nd_smooth(:,3)    - C_nd_smooth(:,4))    / (dx_nd^2);
            d2Cdxx_nd(:,nx) = (2*C_nd_smooth(:,nx) - 5*C_nd_smooth(:,nx-1) + 4*C_nd_smooth(:,nx-2) - C_nd_smooth(:,nx-3)) / (dx_nd^2);
        end
    end

    % ∂²C*/∂y*²: second derivatives (central difference)
    for j = 2:ny-1
        d2Cdyy_nd(j,:) = (C_nd_smooth(j+1,:) - 2*C_nd_smooth(j,:) + C_nd_smooth(j-1,:)) / (dy_nd^2);
    end
    d2Cdyy_nd(1,:)  = (2*C_nd_smooth(1,:)  - 5*C_nd_smooth(2,:)    + 4*C_nd_smooth(3,:)    - C_nd_smooth(4,:))    / (dy_nd^2);
    d2Cdyy_nd(ny,:) = (2*C_nd_smooth(ny,:) - 5*C_nd_smooth(ny-1,:) + 4*C_nd_smooth(ny-2,:) - C_nd_smooth(ny-3,:)) / (dy_nd^2);

    dCdx_nd_all{s}   = dCdx_nd;
    dCdy_nd_all{s}   = dCdy_nd;
    d2Cdxx_nd_all{s}  = d2Cdxx_nd;
    d2Cdyy_nd_all{s}  = d2Cdyy_nd;
end
fprintf('✓ All derivatives computed (5-point stencil)\n\n');

%% ========================================================================
%% SECTION 7: INTEGRAL BALANCE & GAUSSIAN-SHAPED SOURCE TERMS
%% ========================================================================
fprintf('\n=== COMPUTING INTEGRAL BALANCE & SOURCE TERMS (Non-dimensional) ===\n');

% Initialize storage
Ix_all        = cell(1,3);  % integral C dx
dIx_dy_all    = cell(1,3);  % d/dy (integral C dx)
S_int_x_all   = cell(1,3);  % integral S dx from balance
S_gauss_all   = cell(1,3);  % Gaussian-shaped source distribution

for s = 1:3
    fprintf('Processing species: %s...\n', species_names{s});
    C_nd     = C_nd_all{s};        % Ny x Nx
    dCdx_nd  = dCdx_nd_all{s};     % Ny x Nx
    dCdy_nd  = dCdy_nd_all{s};     % Ny x Nx
    d2Cdxx_nd = d2Cdxx_nd_all{s};  % Ny x Nx
    d2Cdyy_nd = d2Cdyy_nd_all{s};  % Ny x Nx

    % === Step 1: Integrate concentration over x ===
    Ix = trapz(x_nd, C_nd, 2);   % Ny x 1
    Ix_all{s} = Ix;

    % === Step 2: Differentiate with respect to y (after smoothing) ===
    Ix_smooth = smooth(Ix, 5); % Smooth before differentiating
    dIx_dy = gradient(Ix_smooth, dy_nd);  % Ny x 1
    dIx_dy_all{s} = dIx_dy;

    % === Step 3: integral of S over x from the width-integrated balance ===
    % Diffusive flux D dC/dx through the two side boundaries
    if nx > 1
        flux_left  = dCdx_nd(:,1);   % at x = left boundary
        flux_right = dCdx_nd(:,end); % at x = right boundary
        boundary_term = flux_right - flux_left; % Ny x 1
    else
        boundary_term = zeros(ny, 1);
    end

    % A uniform D* = 1 (1e-3 mm^2/s) is used for the boundary flux in this step
    D_avg = 1.0;
    diffusion_term = D_avg * boundary_term;

    S_int_x = v_advection_nd .* dIx_dy - diffusion_term;  % Ny x 1
    S_int_x_all{s} = S_int_x;

    % === Step 4: distribute the integral over x with a Gaussian ===
    S_gauss = zeros(ny,nx);
    for j = 1:ny
        % Avoid division by zero
        if abs(S_int_x(j)) > 1e-12
            % Normalized Gaussian along x with proper normalization
            x_range = max(x_nd) - min(x_nd);
            sigma_eff = min(sig_nd, x_range/6); % Limit sigma to reasonable range
            G_norm = 1/(sqrt(2*pi)*sigma_eff);
            G = G_norm * exp(-((x_nd - x_center_S_nd(j)).^2) / (2*sigma_eff^2));
            % Normalize to integrate to 1
            if trapz(x_nd, G) > 1e-12
                G = G / trapz(x_nd, G);
            end
            % Scale by integral S dx at this y
            S_gauss(j,:) = S_int_x(j) .* G;
        end
    end
    S_gauss_all{s} = S_gauss;
end

fprintf('\nSource term calculation completed for all species.\n');
fprintf('  S1 (H2SO4)  range: [%.3e, %.3e]\n', min(S_gauss_all{1}(:)), max(S_gauss_all{1}(:)));
fprintf('  S2 (NaOH)   range: [%.3e, %.3e]\n', min(S_gauss_all{2}(:)), max(S_gauss_all{2}(:)));
fprintf('  S3 (Na2SO4) range: [%.3e, %.3e]\n\n', min(S_gauss_all{3}(:)), max(S_gauss_all{3}(:)));

%% ========================================================================
%% SECTION 8: DIFFUSION ESTIMATION
%% ========================================================================
fprintf('\n=== STEP 2: DIRECT DIFFUSION ESTIMATION ===\n');

D_refined = cell(1,3);

for s = 1:3
    fprintf('Diffusion estimation for %s...\n', species_names{s});

    % Use non-dimensionalized variables
    dCdy_nd = dCdy_nd_all{s};
    dCdx_nd = dCdx_nd_all{s};
    d2Cdxx_nd = d2Cdxx_nd_all{s};
    d2Cdyy_nd = d2Cdyy_nd_all{s};
    S_nd = S_gauss_all{s};

    % Initial estimate using steady-state assumption
    % 0 = D* d2C*/dx*2 - v* dC*/dy* + S*   (plus d2C*/dy*2 if use_full_laplacian)
    if use_full_laplacian, laplacian = d2Cdxx_nd + d2Cdyy_nd; else, laplacian = d2Cdxx_nd; end  % Eq. (2): d2C/dx2 only

    % Initialize D with reasonable bounds
    D_nd = zeros(ny, nx);

    % Only compute where laplacian is significant
    valid_mask = abs(laplacian) > lap_thr;

    if any(valid_mask(:))
        % Compute D where valid
        numerator = v_advection_nd * dCdy_nd - S_nd;
        D_nd(valid_mask) = numerator(valid_mask) ./ laplacian(valid_mask);

        % Apply physical constraints
        D_nd(D_nd < 0) = NaN;      % No negative diffusion
        D_nd(D_nd > D_upper_nd) = NaN;   % Upper bound in non-dimensional units
        D_nd(~valid_mask) = NaN;   % Set invalid regions to NaN

        % Iterative refinement including gradient terms
        for iter = 1:3
            % Smooth current estimate
            D_smooth = D_nd;
            valid_idx = ~isnan(D_nd);
            if sum(valid_idx(:)) > 10
                % Use inpainting to fill NaN values for gradient computation
                % 'EndValues','nearest': no extrapolation beyond the valid points
                D_filled = fillmissing(D_nd, 'linear', 2, 'EndValues', 'nearest');
                D_filled = fillmissing(D_filled, 'linear', 1, 'EndValues', 'nearest');
                D_smooth = imgaussfilt(D_filled, 1.0);
                D_smooth(~valid_idx) = NaN; % Restore NaN where originally invalid
            end

            % Compute dD/dx with finite differences
            dDdx_nd = zeros(ny, nx);
            if nx > 2
                for i = 2:nx-1
                    dDdx_nd(:,i) = (D_smooth(:,i+1) - D_smooth(:,i-1)) / (2*dx_nd);
                end
            end

            % Update D using full PDE
            % 0 = D*Laplacian(C) + grad(D).grad(C) - v*dC/dy* + S*
            grad_D_dot_grad_C = dDdx_nd .* dCdx_nd; % x-component only, as in Eq. (2)
            numerator_refined = v_advection_nd * dCdy_nd - grad_D_dot_grad_C - S_nd;

            D_new = numerator_refined ./ laplacian;

            % Apply constraints
            D_new(D_new < 0) = NaN;
            D_new(D_new > D_upper_nd) = NaN;
            D_new(abs(laplacian) <= lap_thr) = NaN;

            % Update with relaxation
            alpha = 0.3; % Relaxation parameter
            valid_new = ~isnan(D_new);
            valid_old = ~isnan(D_nd);
            update_mask = valid_new & valid_old;

            if any(update_mask(:))
                D_nd(update_mask) = (1-alpha) * D_nd(update_mask) + alpha * D_new(update_mask);
            end
            D_nd(valid_new & ~valid_old) = D_new(valid_new & ~valid_old);
        end

        % Final smoothing while preserving boundaries
        if sum(~isnan(D_nd(:))) > 20
            D_final = D_nd;
            valid_mask_final = ~isnan(D_nd);
            % 'EndValues','nearest': linear extrapolation would create large negative values
            % outside the valid region, which the Gaussian filter would mix into valid points.
            D_filled_final = fillmissing(D_nd, 'linear', 2, 'EndValues', 'nearest');
            D_filled_final = fillmissing(D_filled_final, 'linear', 1, 'EndValues', 'nearest');
            D_smoothed = imgaussfilt(D_filled_final, 0.8);
            D_final(valid_mask_final) = D_smoothed(valid_mask_final);
            D_nd = D_final;
        end
        % re-apply the physical constraints after the smoothing (0 <= D* < D_upper_nd)
        n_out = sum(D_nd(:) < 0 | D_nd(:) > D_upper_nd);
        D_nd(D_nd < 0 | D_nd > D_upper_nd) = NaN;
        fprintf('    Points outside 0 <= D* < %g after smoothing (set to NaN): %d\n', D_upper_nd, n_out);
    else
        fprintf('  Warning: No valid regions for diffusion estimation\n');
        D_nd(:) = NaN;
    end

    D_refined{s} = D_nd*ndD;   % physical units, mm^2/s

    % Summary of the estimate
    valid_D = D_nd(~isnan(D_nd));
    if ~isempty(valid_D)
        fprintf('    Final D_nd range: [%.6f, %.6f]\n', min(valid_D), max(valid_D));
        fprintf('    Final D_nd mean+/-std: %.6f +/- %.6f\n', mean(valid_D), std(valid_D));
        fprintf('    D_phys mean: %.4e mm2/s  (lit: %.4e mm2/s)\n', ...
                mean(valid_D)*ndD, D_lit_array(s));
    else
        fprintf('    Final D_nd range: [NaN, NaN]\n');
    end

    fprintf('    Valid D points: %d/%d (%.1f%%)\n', sum(~isnan(D_nd(:))), numel(D_nd), ...
            100*sum(~isnan(D_nd(:)))/numel(D_nd));
end

%% Save diffusion coefficients into separate variables
D1_refined = D_refined{1}; % H2SO4 diffusion coefficient [mm^2/s]
D2_refined = D_refined{2}; % NaOH  diffusion coefficient [mm^2/s]
D3_refined = D_refined{3}; % Na2SO4 diffusion coefficient [mm^2/s]

fprintf('\nDiffusion estimation completed.\n\n');

%% ========================================================================
%% SECTION 9: SAVE INTERMEDIATE RESULTS
%% ========================================================================
fprintf('=== SAVING INTERMEDIATE RESULTS ===\n');

x_phys = x_phys_mm;
y_phys = y_phys_mm;
X = X_phys;  Y = Y_phys;
C1 = C1_nd;  C2 = C2_nd;  C3 = C3_nd;
S_phys_all = cellfun(@(s) s*ndS, S_gauss_all, 'UniformOutput', false);
S1 = S_phys_all{1};  S2 = S_phys_all{2};  S3 = S_phys_all{3};   % mol/(L s)
dx = dx_phys;  dy = dy_phys;

save('step1_intermediate.mat', ...
     'x_phys', 'y_phys', 'X', 'Y', 'X_phys', 'Y_phys', ...
     'X_nd', 'Y_nd', 'x_nd', 'y_nd', ...
     'x_phys_mm', 'y_phys_mm', ...
     'C1_obs', 'C2_obs', 'C3_obs', ...
     'C1', 'C2', 'C3', 'C_nd_all', 'C_obs_all', ...
     'D1_refined', 'D2_refined', 'D3_refined', 'D_refined', ...
     'S1', 'S2', 'S3', 'S_gauss_all', ...
     'dCdx_nd_all', 'dCdy_nd_all', 'd2Cdxx_nd_all', 'd2Cdyy_nd_all', ...
     'dx', 'dy', 'dx_nd', 'dy_nd', 'dx_phys', 'dy_phys', ...
     'x_min', 'x_max', 'y_max', 'nx', 'ny', ...
     'ndV', 'ndL', 'ndD', 'ndS', 'ndC', ...
     'v_advection_nd', 'v_dimensional', 'sig_nd', ...
     'species_names', 'species_tags', ...
     'D_lit', 'D_lit_array', ...
     'mu1_interp', 'sigma1_interp', 'mu2_interp', 'sigma2_interp', ...
     'A3_interp', 'sigma3_interp', 'condition', 'T_C', 'pixel_mm', 'mu_row');

fprintf('✓ Results saved to: step1_intermediate.mat\n');
fprintf('  → Continuing with the analysis and figures...\n\n');

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                   ESTIMATION COMPLETE                          ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% The workspace is not cleared: the analysis below uses the variables in memory.

%% ========================================================================
%% ANALYSIS, VALIDATION AND FIGURES
%% ========================================================================
%% Mask-threshold sensitivity, validity-mask statistics, step-discontinuity check,
%% cross-validation, comparison with literature values, and figures.
%% ========================================================================

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  STEP 1: ANALYSIS, VALIDATION AND FIGURES                      ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('  Grid: %d × %d, Species: %s, %s, %s\n\n', ...
        ny, nx, species_names{1}, species_names{2}, species_names{3});


%% ========================================================================
%% SECTION 1b: SOURCE TERMS IN PHYSICAL UNITS
%% ========================================================================
% Source terms in physical units (mol/(L s)) for plots and saving;
% S_gauss_all stays non-dimensional for use with the non-dimensional derivatives.
S_phys_all = cellfun(@(s) s*ndS, S_gauss_all, 'UniformOutput', false);
fprintf('=== SOURCE TERMS (physical units, ndS = %.3g mol/(L s)) ===\n', ndS);
for s = 1:3
    fprintf('    %-7s S range: [%.3e, %.3e] mol/(L s)\n', species_tags{s}, min(S_phys_all{s}(:)), max(S_phys_all{s}(:)));
end
fprintf('\n');

%% ========================================================================
%% SECTION 2: SENSITIVITY ANALYSIS FOR MASK THRESHOLDS
%% ========================================================================
fprintf('=== SENSITIVITY ANALYSIS FOR MASK THRESHOLDS ===\n');

laplacian_thresholds     = [0.1, 0.5, 1, 5, 10] * lap_thr;    % around the value used (third entry)
concentration_thresholds = [0.005, 0.01, 0.02, 0.05, 0.10];
gradient_thresholds      = [0.1, 0.5, 1, 5, 10] * grad_thr;

sensitivity_results = struct();

for s = 1:3
    C_nd     = C_nd_all{s};
    dCdx_nd  = dCdx_nd_all{s};
    if use_full_laplacian, laplacian = d2Cdxx_nd_all{s} + d2Cdyy_nd_all{s}; else, laplacian = d2Cdxx_nd_all{s}; end

    cov_lap  = zeros(size(laplacian_thresholds));
    cov_conc = zeros(size(concentration_thresholds));
    cov_grad = zeros(size(gradient_thresholds));

    for k = 1:length(laplacian_thresholds)
        m = abs(laplacian) > laplacian_thresholds(k);
        cov_lap(k) = 100 * sum(m(:)) / numel(m);
    end
    for k = 1:length(concentration_thresholds)
        m = abs(C_nd) > concentration_thresholds(k);
        cov_conc(k) = 100 * sum(m(:)) / numel(m);
    end
    for k = 1:length(gradient_thresholds)
        m = abs(dCdx_nd) > gradient_thresholds(k);
        cov_grad(k) = 100 * sum(m(:)) / numel(m);
    end

    sensitivity_results(s).laplacian     = struct('thresholds', laplacian_thresholds,     'coverage', cov_lap);
    sensitivity_results(s).concentration = struct('thresholds', concentration_thresholds, 'coverage', cov_conc);
    sensitivity_results(s).gradient      = struct('thresholds', gradient_thresholds,      'coverage', cov_grad);

    fprintf('%s:\n', species_names{s});
    fprintf('  Laplacian     threshold %.1e → %.1f%% coverage\n', lap_thr,  cov_lap(3));
    fprintf('  Concentration threshold %.3f → %.1f%% coverage\n', 0.01,  cov_conc(3));
    fprintf('  Gradient      threshold %.1e → %.1f%% coverage\n', grad_thr, cov_grad(3));
end
fprintf('\n');

optimal_thresholds = struct('laplacian', lap_thr, 'concentration', 0.01, 'gradient', grad_thr);

fprintf('Selected thresholds:\n');
fprintf('  |∇²C|     > %.1e\n', optimal_thresholds.laplacian);
fprintf('  |C|       > %.2f\n',  optimal_thresholds.concentration);
fprintf('  |∂C/∂x|   > %.1e\n\n', optimal_thresholds.gradient);

%% ========================================================================
%% SECTION 3: BUILD VALIDITY MASKS AND PRINT STATISTICS
%% ========================================================================
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║          VALIDITY MASK STATISTICS (ALL CRITERIA)              ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

boundary_margin_x = 5;
boundary_margin_y = 10;

valid_masks = cell(1,3);
mask1_all   = cell(1,3);
mask2_all   = cell(1,3);
mask3_all   = cell(1,3);
mask4_ref   = true(ny, nx);
mask4_ref(:, 1:boundary_margin_x)         = false;
mask4_ref(:, end-boundary_margin_x+1:end) = false;
mask4_ref(1:boundary_margin_y, :)         = false;
mask4_ref(end-boundary_margin_y+1:end, :) = false;

n_total = ny * nx;

for s = 1:3
    C_nd      = C_nd_all{s};
    dCdx_nd   = dCdx_nd_all{s};
    if use_full_laplacian, laplacian = d2Cdxx_nd_all{s} + d2Cdyy_nd_all{s}; else, laplacian = d2Cdxx_nd_all{s}; end

    m1 = abs(laplacian) > optimal_thresholds.laplacian;
    m2 = abs(C_nd)      > optimal_thresholds.concentration;
    m3 = abs(dCdx_nd)   > optimal_thresholds.gradient;
    m4 = mask4_ref;
    m_combined = m1 & m2 & m3 & m4;

    mask1_all{s}   = m1;
    mask2_all{s}   = m2;
    mask3_all{s}   = m3;
    valid_masks{s} = m_combined;

    fprintf('%s:\n', species_names{s});
    fprintf('  Criterion 1 (|∇²C|  > %.1e):   %d/%d valid (%.1f%%)\n', lap_thr, ...
            sum(m1(:)), n_total, 100*sum(m1(:))/n_total);
    fprintf('  Criterion 2 (|C|    > 0.01):       %d/%d valid (%.1f%%)\n', ...
            sum(m2(:)), n_total, 100*sum(m2(:))/n_total);
    fprintf('  Criterion 3 (|∂C/∂x|> %.1e):   %d/%d valid (%.1f%%)\n', grad_thr, ...
            sum(m3(:)), n_total, 100*sum(m3(:))/n_total);
    fprintf('  Criterion 4 (boundary exclusion):  %d/%d valid (%.1f%%)\n', ...
            sum(m4(:)), n_total, 100*sum(m4(:))/n_total);
    fprintf('  → COMBINED (all criteria):          %d/%d valid (%.1f%%)\n\n', ...
            sum(m_combined(:)), n_total, 100*sum(m_combined(:))/n_total);
end

%% ========================================================================
%% SECTION 4: STEP DISCONTINUITY CHECK
%% ========================================================================
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║     STEP DISCONTINUITY CHECK (at y ≈ 4 mm)                    ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

[~, y_4mm_idx] = min(abs(y_phys_mm - 4.0));
fprintf('Checking for discontinuities at y ≈ %.2f mm (index %d)\n\n', ...
        y_phys_mm(y_4mm_idx), y_4mm_idx);

for s = 1:3
    D_current    = D_refined{s};
    D_centerline = D_current(:, round(nx/2));
    dD_dy        = gradient(D_centerline, dy_phys);

    window       = max(1, y_4mm_idx-5):min(ny, y_4mm_idx+5);
    max_gradient = max(abs(dD_dy(window)));
    mean_gradient = mean(abs(dD_dy));

    fprintf('%s diffusion coefficient:\n', species_names{s});
    fprintf('  Max |dD/dy| near y=4mm: %.6e\n', max_gradient);
    fprintf('  Mean |dD/dy| overall:   %.6e\n', mean_gradient);
    fprintf('  Ratio: %.2f\n', max_gradient / (mean_gradient + 1e-12));
    if max_gradient > 5 * mean_gradient
        fprintf('  ⚠ DISCONTINUITY DETECTED!\n');
    else
        fprintf('  ✓ Smooth (no discontinuity)\n');
    end
    fprintf('\n');
end

for s = 1:3
    S_current    = S_gauss_all{s};
    S_centerline = S_current(:, round(nx/2));
    dS_dy        = gradient(S_centerline, dy_phys);

    window       = max(1, y_4mm_idx-5):min(ny, y_4mm_idx+5);
    max_gradient = max(abs(dS_dy(window)));
    mean_gradient = mean(abs(dS_dy));

    fprintf('%s source term:\n', species_names{s});
    fprintf('  Max |dS/dy| near y=4mm: %.6e\n', max_gradient);
    fprintf('  Mean |dS/dy| overall:   %.6e\n', mean_gradient);
    fprintf('  Ratio: %.2f\n', max_gradient / (mean_gradient + 1e-12));
    if max_gradient > 5 * mean_gradient
        fprintf('  ⚠ DISCONTINUITY DETECTED!\n');
    else
        fprintf('  ✓ Smooth (no discontinuity)\n');
    end
    fprintf('\n');
end

%% ========================================================================
%% SECTION 5: CROSS-VALIDATION
%% ========================================================================
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                   CROSS-VALIDATION                            ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

K_folds    = 5;
lambda_base = 0.01;
min_points  = 10;
cv_results  = cell(1,3);

for s = 1:3
    fprintf('%s:\n', species_names{s});

    C_nd      = C_nd_all{s};
    dCdy_nd   = dCdy_nd_all{s};
    if use_full_laplacian, laplacian = d2Cdxx_nd_all{s} + d2Cdyy_nd_all{s}; else, laplacian = d2Cdxx_nd_all{s}; end
    S_nd      = S_gauss_all{s};

    valid_mask         = valid_masks{s};
    [valid_j, valid_i] = find(valid_mask);
    n_valid            = length(valid_j);

    if n_valid < K_folds * min_points
        fprintf('  Insufficient points for cross-validation\n\n');
        continue;
    end

    rng(42);
    fold_assignment = mod(randperm(n_valid), K_folds) + 1;
    cv_errors       = zeros(K_folds, 1);
    cv_D_estimates  = zeros(K_folds, 1);

    for k = 1:K_folds
        train_idx = fold_assignment ~= k;
        test_idx  = fold_assignment == k;

        idx_tr = sub2ind([ny,nx], valid_j(train_idx), valid_i(train_idx));
        A_tr   = laplacian(idx_tr);
        b_tr   = v_advection_nd * dCdy_nd(idx_tr) - S_nd(idx_tr);

        w  = sqrt(abs(A_tr) / (sum(abs(A_tr)) + 1e-10));
        Aw = A_tr .* w;  bw = b_tr .* w;
        D_k = (Aw' * bw) / (Aw' * Aw + lambda_base);
        cv_D_estimates(k) = D_k;

        idx_te    = sub2ind([ny,nx], valid_j(test_idx), valid_i(test_idx));
        A_te      = laplacian(idx_te);
        b_te      = v_advection_nd * dCdy_nd(idx_te) - S_nd(idx_te);
        residuals = b_te - D_k * A_te;
        cv_errors(k) = sqrt(mean(residuals.^2));
    end

    cv_results{s} = struct( ...
        'D_estimates', cv_D_estimates, ...
        'errors',      cv_errors, ...
        'mean_D',      mean(cv_D_estimates), ...
        'std_D',       std(cv_D_estimates), ...
        'mean_RMSE',   mean(cv_errors));

    fprintf('  Cross-validated D: %.4e ± %.4e\n', ...
            mean(cv_D_estimates), std(cv_D_estimates));
    fprintf('  Mean RMSE:         %.4e\n', mean(cv_errors));
    fprintf('  CV coeff. of variation: %.2f%%\n\n', ...
            100*std(cv_D_estimates)/mean(cv_D_estimates));
end

%% ========================================================================
%% SECTION 6: COMPARISON WITH LITERATURE VALUES
%% ========================================================================
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║              COMPARISON WITH LITERATURE                       ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

convergence_info = cell(1,3);

fprintf('%-10s | %12s | %12s | %8s | %10s\n', ...
        'Species', 'Estimated', 'Literature', 'Ratio', 'Status');
fprintf('%s\n', repmat('-', 1, 60));

for s = 1:3
    valid_D = D_refined{s}(~isnan(D_refined{s}));
    if ~isempty(valid_D)
        D_est   = mean(valid_D);
        D_lit_v = D_lit_array(s);
        ratio   = D_est / D_lit_v;
        if     ratio > 0.5 && ratio < 2.0, status = '✓ Good';
        elseif ratio > 0.2 && ratio < 5.0, status = '~ Fair';
        else,                               status = '✗ Check';
        end
        fprintf('%-10s | %10.2e   | %10.2e   | %8.2f | %10s\n', ...
                species_tags{s}, D_est, D_lit_v, ratio, status);

        convergence_info{s} = struct( ...
            'mean',      mean(valid_D), ...
            'std',       std(valid_D), ...
            'literature', D_lit_array(s), ...
            'ratio',     mean(valid_D)/D_lit_array(s), ...
            'coverage',  100*length(valid_D)/numel(D_refined{s}));
    else
        fprintf('%-10s | %12s | %10.2e   | %8s | %10s\n', ...
                species_tags{s}, 'N/A', D_lit_array(s), 'N/A', 'No data');
    end
end
fprintf('\n');

%% ========================================================================
%% SECTION 7: FIGURES
%% ------------------------------------------------------------------------
%% ========================================================================
fprintf('═══════════════════════════════════════════════════════════\n');
fprintf('=== CREATING VISUALIZATIONS ===\n');
fprintf('═══════════════════════════════════════════════════════════\n\n');

make_mask_rgb = @(m) cat(3, double(~m), double(m), double(m));

% Standard label strings
LX = 'x (mm)';
LY = 'y (mm)';
LC = 'C (mol/L)';
LD = 'D (mm^2/s)';
LS = 'S (mol/(L\cdot s))';

%% ── Figure 1: Profile Parameter Evolution ──────────────────────────────
figure('Name', 'Profile Parameter Evolution', 'Position', [50, 50, 1400, 800]);

mu1_profile    = mu1_interp(:, round(nx/2));
sigma1_profile = sigma1_interp(:, round(nx/2));
mu2_profile    = mu2_interp(:, round(nx/2));
sigma2_profile = sigma2_interp(:, round(nx/2));
A3_profile     = A3_interp(:, round(nx/2));
sigma3_profile = sigma3_interp(:, round(nx/2));
dmu1_dy = gradient(mu1_profile, dy_nd);
dmu2_dy = gradient(mu2_profile, dy_nd);

subplot(2,3,1);
plot(y_phys_mm, mu1_profile, 'b-', 'LineWidth', 2);
xlabel(LY); ylabel('\mu_1  (center, non-dim)'); title('H_2SO_4 center'); grid on;

subplot(2,3,2);
plot(y_phys_mm, sigma1_profile, 'b-', 'LineWidth', 2);
xlabel(LY); ylabel('\sigma_1  (width, non-dim)'); title('H_2SO_4 width'); grid on;

subplot(2,3,3);
plot(y_phys_mm, A3_profile, 'g-', 'LineWidth', 2); hold on;
plot(y_phys_mm, sigma3_profile, 'g--', 'LineWidth', 2); hold off;
xlabel(LY); ylabel('amplitude / width');
title('Na_2SO_4 Gaussian parameters');
legend({'A_3 (amplitude)', '\sigma_3 (width)'}, 'Location', 'best'); grid on;

subplot(2,3,4);
plot(y_phys_mm, mu2_profile, 'r-', 'LineWidth', 2);
xlabel(LY); ylabel('\mu_2  (center, non-dim)'); title('NaOH center'); grid on;

subplot(2,3,5);
plot(y_phys_mm, sigma2_profile, 'r-', 'LineWidth', 2);
xlabel(LY); ylabel('\sigma_2  (width, non-dim)'); title('NaOH width'); grid on;

subplot(2,3,6);
plot(y_phys_mm, dmu1_dy, 'b-', 'LineWidth', 2); hold on;
plot(y_phys_mm, dmu2_dy, 'r-', 'LineWidth', 2); hold off;
xlabel(LY); ylabel('d\mu/dy'); title('Center-position gradient');
legend({'H_2SO_4', 'NaOH'}, 'Location', 'best'); grid on;

%% ── Figure 2: Threshold Sensitivity ────────────────────────────────────
figure('Name', 'Threshold Sensitivity', 'Position', [100, 100, 1400, 400]);
for s = 1:3
    subplot(1,3,s); hold on;
    plot(log10(sensitivity_results(s).laplacian.thresholds), ...
         sensitivity_results(s).laplacian.coverage,  'b-o', 'LineWidth', 2);
    plot(log10(sensitivity_results(s).gradient.thresholds), ...
         sensitivity_results(s).gradient.coverage,   'r-s', 'LineWidth', 2);
    plot(log10(sensitivity_results(s).concentration.thresholds), ...
         sensitivity_results(s).concentration.coverage, 'g-^', 'LineWidth', 2);
    xline(log10(optimal_thresholds.laplacian), 'k--', 'LineWidth', 1.5);
    xlabel('log_{10}(threshold)'); ylabel('valid coverage (%)');
    title(species_names{s});
    if s == 1
        legend({'|\nabla^2C|', '|\partialC/\partialx|', '|C|', 'selected'}, ...
               'Location', 'best');
    end
    grid on; hold off;
end
sgtitle('Mask-Threshold Sensitivity', 'FontWeight', 'bold');

%% ── Figure 3: Validity Masks (top view) ────────────────────────────────
figure('Name', 'Validity Masks (green = valid)', 'Position', [200, 100, 1400, 450]);
for s = 1:3
    subplot(1,3,s);
    imagesc(x_phys_mm, y_phys_mm, make_mask_rgb(valid_masks{s}));
    axis xy; hold on;
    xL = x_phys_mm(boundary_margin_x);
    xR = x_phys_mm(end - boundary_margin_x + 1);
    yB = y_phys_mm(boundary_margin_y);
    yT = y_phys_mm(end - boundary_margin_y + 1);
    plot([xL xL], [y_phys_mm(1) y_phys_mm(end)], 'r--', 'LineWidth', 1.5);
    plot([xR xR], [y_phys_mm(1) y_phys_mm(end)], 'r--', 'LineWidth', 1.5);
    plot([x_phys_mm(1) x_phys_mm(end)], [yB yB], 'r--', 'LineWidth', 1.5);
    plot([x_phys_mm(1) x_phys_mm(end)], [yT yT], 'r--', 'LineWidth', 1.5);
    hold off;
    xlabel(LX); ylabel(LY);
    cov_v = 100 * sum(valid_masks{s}(:)) / n_total;
    title(sprintf('%s  (%.0f%% valid)', species_names{s}, cov_v));
end
sgtitle('Validity Masks  —  red dashed = boundary exclusion', 'FontWeight', 'bold');

%% ── Figure 4: Mask Criteria Breakdown ──────────────────────────────────
criteria_labels = { '|\nabla^2C| threshold', '|C| threshold', ...
                    '|\partialC/\partialx| threshold', 'boundary exclusion', ...
                    'combined valid' };

figure('Name', 'Mask Criteria Breakdown (green = valid)', 'Position', [300, 150, 1400, 900]);
for s = 1:3
    masks_array = {mask1_all{s}, mask2_all{s}, mask3_all{s}, mask4_ref, valid_masks{s}};
    for m = 1:5
        subplot(3, 5, (s-1)*5 + m);
        imagesc(x_phys_mm, y_phys_mm, make_mask_rgb(masks_array{m}));
        axis xy;
        if s == 1, title(criteria_labels{m}, 'Interpreter', 'tex', 'FontSize', 9); end
        if m == 1, ylabel(sprintf('%s\n%s', species_names{s}, LY));
        else,      ylabel(LY); end
        if s == 3, xlabel(LX); end
        cov_m = 100 * sum(masks_array{m}(:)) / n_total;
        subtitle(sprintf('%.0f%%', cov_m), 'FontSize', 8);
    end
end
sgtitle('Validity-Mask Criteria (green = valid, per-panel % shown)', 'FontWeight', 'bold');

%% ── Figure 5: Concentration Profiles ───────────────────────────────────
[~, y_4mm_idx] = min(abs(y_phys_mm - 4.0));

figure('Name', 'Concentration Profiles (ERFC / Gaussian)', 'Position', [100, 50, 1400, 900]);
for s = 1:3
    subplot(2,3,s);
    surf(X_phys, Y_phys, C_nd_all{s});
    shading interp; colormap parula; c = colorbar; ylabel(c, LC);
    xlabel(LX); ylabel(LY); zlabel(LC);
    title(species_names{s}); view(45, 30);

    subplot(2,3,s+3);
    plot(y_phys_mm, C_nd_all{s}(:, round(nx/2)), 'b-', 'LineWidth', 2);
    xlabel(LY); ylabel(LC);
    title(sprintf('%s  centerline', species_names{s})); grid on;
    hold on; plot([4 4], ylim, 'r--', 'LineWidth', 1.5); hold off;
end
sgtitle('Modelled Concentration Fields and Centerlines', 'FontWeight', 'bold');

%% ── Figure 6: Diffusion Coefficients (mask overlay) ────────────────────
figure('Name', 'Diffusion Coefficients (red contour = masked)', ...
       'Position', [250, 250, 1400, 450]);
for s = 1:3
    subplot(1,3,s);
    surf(X_phys, Y_phys, D_refined{s});
    shading interp; hold on;
    contour3(X_phys, Y_phys, double(~valid_masks{s}), [0.5 0.5], ...
             'LineColor', 'r', 'LineWidth', 2);
    hold off;
    colormap parula; c = colorbar; ylabel(c, LD);
    xlabel(LX); ylabel(LY); zlabel(LD);
    if ~isempty(convergence_info{s})
        title(sprintf('%s  (D_{est}/D_{lit} = %.2f)', ...
              species_names{s}, convergence_info{s}.ratio));
    else
        title(species_names{s});
    end
    view(45, 30);
end
sgtitle('Estimated Diffusion Coefficients  —  red contour = masked region', ...
        'FontWeight', 'bold');

%% ── Figure 7: Observed Concentrations (mask overlay) ───────────────────
figure('Name', 'Observed Concentrations (red contour = masked)', ...
       'Position', [100, 50, 1400, 450]);
for s = 1:3
    subplot(1,3,s);
    surf(X_phys, Y_phys, C_obs_all{s});
    shading interp; hold on;
    contour3(X_phys, Y_phys, double(~valid_masks{s}), [0.5 0.5], ...
             'LineColor', 'r', 'LineWidth', 2);
    hold off;
    colormap parula; c = colorbar; ylabel(c, LC);
    xlabel(LX); ylabel(LY); zlabel(LC);
    title(species_names{s}); view(45, 30);
end
sgtitle('Observed Concentration Fields  —  red contour = masked region', ...
        'FontWeight', 'bold');

%% ── Figure 8: Modelled Concentrations (mask overlay) ───────────────────
figure('Name', 'Modelled Concentrations (red contour = masked)', ...
       'Position', [100, 500, 1400, 450]);
for s = 1:3
    subplot(1,3,s);
    surf(X_phys, Y_phys, C_nd_all{s});
    shading interp; hold on;
    contour3(X_phys, Y_phys, double(~valid_masks{s}), [0.5 0.5], ...
             'LineColor', 'r', 'LineWidth', 2);
    hold off;
    colormap parula; c = colorbar; ylabel(c, LC);
    xlabel(LX); ylabel(LY); zlabel(LC);
    title(species_names{s}); view(45, 30);
end
sgtitle('Modelled Concentration Fields  —  red contour = masked region', ...
        'FontWeight', 'bold');

%% ── Figure 9: Diffusion Coefficients (centerline check) ────────────────
figure('Name', 'Diffusion Coefficients + Centerline', ...
       'Position', [200, 200, 1400, 900]);
for s = 1:3
    subplot(2,3,s);
    surf(X_phys, Y_phys, D_refined{s});
    shading interp; colormap parula; c = colorbar; ylabel(c, LD);
    xlabel(LX); ylabel(LY); zlabel(LD);
    title(species_names{s}); view(45, 30);

    subplot(2,3,s+3);
    plot(y_phys_mm, D_refined{s}(:, round(nx/2)), 'b-', 'LineWidth', 2);
    xlabel(LY); ylabel(LD);
    title(sprintf('%s  centerline', species_names{s})); grid on;
    hold on; plot([4 4], ylim, 'r--', 'LineWidth', 1.5); hold off;
end
sgtitle('Estimated Diffusion Fields and Centerlines', 'FontWeight', 'bold');

%% ── Figure 10: Source Terms ────────────────────────────────────────────
figure('Name', 'Source Terms', 'Position', [350, 350, 1400, 450]);
for s = 1:3
    subplot(1,3,s);
    surf(X_phys, Y_phys, S_phys_all{s});
    shading interp; colormap parula; c = colorbar; ylabel(c, LS);
    xlabel(LX); ylabel(LY); zlabel(LS);
    title(species_names{s}); view(45, 30);
end
sgtitle('Reconstructed Source Terms', 'FontWeight', 'bold');

%% ── Figure 11: Cross-Validation ────────────────────────────────────────
figure('Name', 'Cross-Validation', 'Position', [200, 200, 1200, 400]);
for s = 1:3
    subplot(1,3,s);
    if ~isempty(cv_results{s})
        bar(cv_results{s}.D_estimates); hold on;
        yline(cv_results{s}.mean_D, 'b-',  'LineWidth', 2);
        yline(D_lit_array(s),       'r--', 'LineWidth', 2);
        hold off;
        xlabel('fold'); ylabel(LD);
        title(species_names{s});
        if s == 1
            legend({'fold estimate', 'CV mean', 'literature'}, 'Location', 'best');
        end
        grid on;
    end
end
sgtitle('K-Fold Cross-Validated Diffusion Estimates', 'FontWeight', 'bold');

%% ── Figure 12: Estimated vs Literature ─────────────────────────────────
figure('Name', 'Estimated vs Literature', 'Position', [300, 300, 600, 500]);
D_est_means = zeros(1,3); D_est_stds = zeros(1,3);
for s = 1:3
    v = D_refined{s}(~isnan(D_refined{s}));
    if ~isempty(v), D_est_means(s) = mean(v); D_est_stds(s) = std(v); end
end
bw = 0.35;
h1 = bar((1:3) - bw/2, D_est_means, bw); set(h1, 'FaceColor', [0.2 0.6 0.8]); hold on;
h2 = errorbar((1:3) - bw/2, D_est_means, D_est_stds, 'k.'); set(h2, 'LineWidth', 1.5);
h3 = bar((1:3) + bw/2, D_lit_array,  bw); set(h3, 'FaceColor', [0.8 0.4 0.2]); hold off;
set(gca, 'XTick', 1:3, 'XTickLabel', species_tags);
ylabel(LD); title('Estimated vs Literature');
legend([h1, h3], {'estimated (mean \pm std)', 'literature'}, 'Location', 'best');
grid on;

%% ── Figure 13: Observed vs Modelled ────────────────────────────────────
figure('Name', 'Observed vs Modelled', 'Position', [400, 400, 1400, 900]);
for s = 1:3
    subplot(3,2,2*s-1);
    surf(X_phys, Y_phys, C_obs_all{s});
    shading interp; colormap parula; c = colorbar; ylabel(c, LC);
    xlabel(LX); ylabel(LY); zlabel(LC);
    title(sprintf('%s  (observed)', species_names{s})); view(45, 30);
    xlim([x_min, x_max]); ylim([0, y_max]);

    subplot(3,2,2*s);
    surf(X_phys, Y_phys, C_nd_all{s});
    shading interp; colormap parula; c = colorbar; ylabel(c, LC);
    xlabel(LX); ylabel(LY); zlabel(LC);
    title(sprintf('%s  (modelled)', species_names{s})); view(45, 30);
    xlim([x_min, x_max]); ylim([0, y_max]);
end
sgtitle('Observed vs Modelled Concentration Fields', 'FontWeight', 'bold');

%% ── Figure 14: Cross-Sectional Profiles ────────────────────────────────
y_sections_mm = [1, 4, 8];
y_section_indices = zeros(size(y_sections_mm));
for i = 1:numel(y_sections_mm)
    [~, y_section_indices(i)] = min(abs(y_phys_mm - y_sections_mm(i)));
end
figure('Name', 'Transverse Profiles by Depth', 'Position', [500, 500, 1200, 800]);
for s = 1:3
    subplot(3,1,s); hold on;
    for i = 1:numel(y_sections_mm)
        idx = y_section_indices(i);
        plot(x_phys_mm, C_nd_all{s}(idx,:), '-', 'LineWidth', 2, ...
             'DisplayName', sprintf('modelled, y = %.0f mm', y_sections_mm(i)));
        plot(x_phys_mm, C_obs_all{s}(idx,:), '--', 'LineWidth', 1.5, ...
             'DisplayName', sprintf('observed, y = %.0f mm', y_sections_mm(i)));
    end
    xlabel(LX); ylabel(LC); title(species_names{s});
    legend('Location', 'eastoutside'); grid on;
    xlim([x_min, x_max]); hold off;
end
sgtitle('Transverse Concentration Profiles at y = 1, 4, 8 mm', 'FontWeight', 'bold');

fprintf('✓ All standard figures created\n\n');

%% ========================================================================
%% SECTION 7b: 2D LINE PROFILE PLOTS
%% ========================================================================
fprintf('═══════════════════════════════════════════════════════════\n');
fprintf('=== CREATING 2D LINE PROFILE PLOTS ===\n');
fprintf('═══════════════════════════════════════════════════════════\n\n');

create_line_profiles_all( ...
    X_phys, Y_phys, x_phys_mm, y_phys_mm, ...
    C_nd_all, C_obs_all, D_refined, S_phys_all, ...
    species_names, species_tags, D_lit_array, valid_masks);

fprintf('✓ All line profile plots created and saved\n\n');
%% ========================================================================
%% SECTION 7c: 2D MAPS OF THE MEASURED CONCENTRATION FIELDS
%% ------------------------------------------------------------------------
%% imagesc plots C_show{s}', so the horizontal axis is y (streamwise) and the
%% vertical axis is x (transverse).
%% ========================================================================
fprintf('═══════════════════════════════════════════════════════════\n');
fprintf('=== CREATING 2D MAPS OF THE MEASURED FIELDS ===\n');
fprintf('═══════════════════════════════════════════════════════════\n\n');

% For display only (Fig. 2(b) of the article), the pixel maps are smoothed with a 3 x 3 median
% filter before the same centring and interpolation as C_obs; the estimation uses C_obs unfiltered.
C_show = cell(1, 3);
for s = 1:3
    Fm = medfilt2(F{s}, [3 3], 'symmetric');
    C_show{s} = zeros(ny, nx);
    for j = 1:ny
        C_show{s}(j,:) = interp1(r_img, Fm(j,:), mu_row(j) + x_phys_mm / pixel_mm, 'linear');
    end
    C_show{s}(C_show{s} < 0) = 0;
end
show_label = 'Observed, 3 x 3 median filter';
C1_med = C_show{1};  C2_med = C_show{2};  C3_med = C_show{3};   % H2SO4, NaOH, Na2SO4 (mol/L)
save('map_concentration_medfilt.mat', 'C1_med', 'C2_med', 'C3_med', 'x_phys_mm', 'y_phys_mm', 'condition');
fprintf('   Saved -> map_concentration_medfilt.mat  (median-filtered maps of Fig. 2(b))\n');

clim_list  = { [0 1.2], [0 2.4], [0 0.6] };
ctick_list = { [0 0.6 1.2], [0 1.2 2.4], [0 0.3 0.6] };
cb_label_list = { 'C_{H_2SO_4} (mol/L)', 'C_{NaOH} (mol/L)', 'C_{Na_2SO_4} (mol/L)' };

y_station_mm = [1.0, 4.0, 8.0];
line_color   = 'w';

x_lo = min(x_phys_mm);  x_hi = max(x_phys_mm);
y_ticks = 0:2:floor(max(y_phys_mm));
x_ticks = round(linspace(x_lo, x_hi, 3), 2);

figT = figure('Name', '2D Concentration Maps', ...
    'Position', [80, 120, 1500, 440], 'Color', 'w');

% Panel geometry (left margin reserved for the inlet arms)
left0 = 0.15; panel_w = 0.245; gap = 0.028;
panel_h = 0.48; panel_b = 0.42; cb_h = 0.045; cb_gap = 0.19;

ax_map = gobjects(1,3);
for s = 1:3
    axL = left0 + (s-1)*(panel_w + gap);
    ax_map(s) = axes('Position', [axL, panel_b, panel_w, panel_h]); %#ok<LAXES>
    imagesc(ax_map(s), y_phys_mm, x_phys_mm, C_show{s}');
    set(ax_map(s), 'YDir', 'normal'); axis(ax_map(s), 'tight');
    colormap(ax_map(s), parula); caxis(ax_map(s), clim_list{s});

    % Title: species name
    title(ax_map(s), species_names{s}, ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'tex');

    set(ax_map(s), 'Box', 'on', 'Layer', 'top', 'LineWidth', 0.8, ...
        'FontSize', 9, 'XTick', y_ticks, 'YTick', x_ticks, ...
        'TickDir', 'out', 'TickLength', [0.015 0]);
    xlabel(ax_map(s), 'y (mm)', 'FontSize', 10);
    if s == 1, ylabel(ax_map(s), 'x (mm)', 'FontSize', 10); end

    hold(ax_map(s), 'on');
    for k = 1:numel(y_station_mm)
        yk = y_station_mm(k);
        plot(ax_map(s), [yk yk], [x_lo x_hi], '--', ...
            'Color', line_color, 'LineWidth', 1.3);
    end
    hold(ax_map(s), 'off');

    % Individual colorbar underneath, labelled per species
    axcb = axes('Position', [axL, panel_b - cb_gap, panel_w, cb_h]); %#ok<LAXES>
    ncol = 256;
    imagesc(axcb, linspace(clim_list{s}(1), clim_list{s}(2), ncol), 1, 1:ncol);
    colormap(axcb, parula);
    set(axcb, 'YDir', 'normal', 'YTick', [], 'XLim', clim_list{s}, ...
        'XTick', ctick_list{s}, 'FontSize', 9);
    axcb.XAxis.TickLength = [0.02 0]; box(axcb, 'on');
    xlabel(axcb, cb_label_list{s}, 'FontSize', 10, 'Interpreter', 'tex');
end

% ------------------------------------------------------------------------
% Coordinate and inlet cues, drawn in a transparent overlay.
% ------------------------------------------------------------------------
axOv = axes('Position', [0 0 1 1], 'Color', 'none'); %#ok<LAXES>
axis(axOv, [0 1 0 1]); axis(axOv, 'off'); hold(axOv, 'on'); set(axOv, 'YDir', 'normal');



hold(axOv, 'off'); uistack(axOv, 'top');

save_figure_safe(figT, 'map_concentration.png', 600);
fprintf('   Saved -> map_concentration.png  (%s field)\n\n', show_label);

%% ========================================================================
%% SECTION 7d: 2D MAPS, MEASURED vs FITTED
%% ------------------------------------------------------------------------
%% One row for the measured maps and one for the fits, with a shared colour scale per species.
%% ========================================================================
fprintf('═══════════════════════════════════════════════════════════\n');
fprintf('=== CREATING 2D MAPS: OBSERVED vs MODELLED ===\n');
fprintf('═══════════════════════════════════════════════════════════\n\n');

figOF = figure('Name', '2D Maps: Observed vs Modelled', ...
               'Position', [60, 80, 1500, 780], 'Color', 'w');

left0 = 0.16; panel_w = 0.235; gap = 0.025; panel_h = 0.30;
row_b = [0.56, 0.10]; cb_h = 0.030; cb_gap = 0.055;

field_rows = {C_obs_all, C_nd_all};
row_labels = {'Observed', 'Modelled'};
ax_map2 = gobjects(2,3);

for r = 1:2
    C_show_row = field_rows{r};
    for s = 1:3
        axL = left0 + (s-1)*(panel_w + gap);
        ax_map2(r,s) = axes('Position', [axL, row_b(r), panel_w, panel_h]); %#ok<LAXES>
        imagesc(ax_map2(r,s), y_phys_mm, x_phys_mm, C_show_row{s}');
        set(ax_map2(r,s), 'YDir', 'normal'); axis(ax_map2(r,s), 'tight');
        colormap(ax_map2(r,s), parula); caxis(ax_map2(r,s), clim_list{s});
        set(ax_map2(r,s), 'XTick', [], 'YTick', [], 'Box', 'on', ...
            'LineWidth', 0.8, 'Layer', 'top');
        if r == 1
            title(ax_map2(r,s), species_names{s}, 'FontSize', 13, ...
                  'FontWeight', 'bold', 'Interpreter', 'tex');
        end
        axcb = axes('Position', [axL, row_b(r) - cb_gap, panel_w, cb_h]); %#ok<LAXES>
        ncol = 256;
        imagesc(axcb, linspace(clim_list{s}(1), clim_list{s}(2), ncol), 1, 1:ncol);
        colormap(axcb, parula);
        set(axcb, 'YDir', 'normal', 'YTick', [], 'XLim', clim_list{s}, ...
            'XTick', ctick_list{s}, 'FontSize', 8);
        axcb.XAxis.TickLength = [0.02 0]; box(axcb, 'on');
        xlabel(axcb, 'C (mol/L)', 'FontSize', 8);
    end
end

axOv2 = axes('Position', [0 0 1 1], 'Color', 'none'); %#ok<LAXES>
axis(axOv2, [0 1 0 1]); axis(axOv2, 'off'); hold(axOv2, 'on'); set(axOv2, 'YDir', 'normal');

hold(axOv2, 'off'); uistack(axOv2, 'top');

save_figure_safe(figOF, 'map_concentration_obs_vs_fit.png', 600);
fprintf('   Saved -> map_concentration_obs_vs_fit.png\n\n');

%% ========================================================================
%% SECTION 7e: 3D SURFACE PLOTS OF D AND S
%% ------------------------------------------------------------------------
%% Row 1: D(x,y) per species, transparent plane = literature value
%% Row 2: S(x,y) per species, transparent plane = S = 0
%% Black lines = sampling stations y = 1, 4, 8 mm (as in the line-profile figures)
%% ========================================================================
fprintf('═══════════════════════════════════════════════════════════\n');
fprintf('=== CREATING 3D PROFILES OF D AND S ===\n');
fprintf('═══════════════════════════════════════════════════════════\n\n');

opts3d = struct();
opts3d.mask     = false;        % true: show only points inside the validity masks
opts3d.robust   = true;         % true: colour/z limits from 1st-99th percentile (hides spikes)
opts3d.stations = [1 4 8];      % y stations (mm) drawn on the surfaces; [] = none
opts3d.view     = [-37.5 28];   % [azimuth elevation]
opts3d.dpi      = 600;

create_3d_profiles_DS(X_phys, Y_phys, D_refined, S_phys_all, ...
    species_names, D_lit_array, valid_masks, opts3d);

fprintf('✓ 3D profiles saved -> profiles_3D_DS.png\n\n');

%% ========================================================================
%% SECTION 8: SAVE FINAL RESULTS
%% ========================================================================
fprintf('=== SAVING FINAL RESULTS ===\n');

x_phys = x_phys_mm;
y_phys = y_phys_mm;
X = X_phys;  Y = Y_phys;
C1 = C_nd_all{1};  C2 = C_nd_all{2};  C3 = C_nd_all{3};
D1_refined = D_refined{1};
D2_refined = D_refined{2};
D3_refined = D_refined{3};
S1 = S_phys_all{1};   % mol/(L s)
S2 = S_phys_all{2};
S3 = S_phys_all{3};
dx = dx_phys;  dy = dy_phys;

S_units = 'mol/(L s)';  D_units = 'mm^2/s';  C_units = 'mol/L';
save('step1_initial_estimates.mat', 'S_units', 'D_units', 'C_units', 'fit_info', 'profile_source', 'x_map_mode', ...
     'x_phys', 'y_phys', 'X', 'Y', ...
     'C1_obs', 'C2_obs', 'C3_obs', ...
     'D1_refined', 'D2_refined', 'D3_refined', ...
     'C1', 'C2', 'C3', ...
     'S1', 'S2', 'S3', ...
     'dx', 'dy', 'x_min', 'x_max', ...
     'ndV', 'ndL', 'ndD', 'ndS', 'ndC', ...
     'species_names', ...
     'mu1_interp', 'sigma1_interp', 'mu2_interp', 'sigma2_interp', ...
     'A3_interp', 'sigma3_interp', ...
     'condition', 'T_C', 'v_dimensional', 'pixel_mm', 'mu_row', 'D_lit_array');

fprintf('✓ Final results saved to: step1_initial_estimates.mat\n\n');

%% ========================================================================
%% FINAL SUMMARY
%% ========================================================================
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                     FINAL SUMMARY                             ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('ESTIMATION QUALITY:\n');
for s = 1:3
    if ~isempty(convergence_info{s})
        fprintf('  %s: D = %.2e mm²/s (%.0f%% of literature, %.1f%% coverage)\n', ...
                species_names{s}, convergence_info{s}.mean, ...
                100*convergence_info{s}.ratio, convergence_info{s}.coverage);
    end
end

fprintf('\nFILES SAVED:\n');
fprintf('  step1_intermediate.mat\n');
fprintf('  step1_initial_estimates.mat\n');
fprintf('  line_profiles_C_transverse.png\n');
fprintf('  line_profiles_C_streamwise.png\n');
fprintf('  line_profiles_D_lines.png\n');
fprintf('  line_profiles_S_lines.png\n');
fprintf('  line_profiles_DS.png\n');
fprintf('  map_concentration.png\n');
fprintf('  map_concentration_obs_vs_fit.png\n');
fprintf('  profiles_3D_DS.png\n');

fprintf('\n✓ Analysis complete! Results are in %s\n', out_dir);
cd(start_dir);
fprintf('Note: Due to ndC = 1, computed values are in physical units directly.\n\n');

fprintf('═══════════════════════════════════════════════════════════\n');

%% ========================================================================
%% LOCAL FUNCTION: LINE PROFILE PLOTS
%% ------------------------------------------------------------------------
%% Local functions must be placed at the end of a script file (MATLAB R2016b or newer).
%%
%% ========================================================================

function create_line_profiles_all(X, Y, x_vec, y_vec, ...
        C_th_all, C_obs_all, D_all, S_all, ...
        species_names, species_tags, D_lit_array, valid_masks)

    if nargin < 12, valid_masks = {}; end %#ok<NASGU>

    [Ny, Nx] = size(X);

    % --- Sampling locations ---------------------------------------------
    y_idx_list   = round(linspace(round(0.08*Ny), round(0.92*Ny), 5));
    x_idx_center = round(Nx/2);

    % --- Style ----------------------------------------------------------
    sp_colors = {[0.85 0.10 0.10], [0.10 0.30 0.85], [0.10 0.65 0.20]};
    sp_marker = {'o', 's', 'd'};
    unit_C    = 'mol/L';
    unit_D    = 'mm^2/s';
    unit_S    = 'mol/(L\cdot s)';

    %% ===================================================================
    %% FIG 1 — Transverse C(x) at multiple y, modelled vs observed
    %% ===================================================================
    fig1 = figure('Name', 'Transverse Concentration Profiles', ...
                  'Position', [50, 50, 1700, 850], 'Color', 'w');
    n_cuts = numel(y_idx_list);
    for sp = 1:3
        for c = 1:n_cuts
            subplot(3, n_cuts, (sp-1)*n_cuts + c);
            hold on; grid on; box on;
            yi = y_idx_list(c);

            plot(x_vec, C_th_all{sp}(yi, :), '-', ...
                'Color', sp_colors{sp}, 'LineWidth', 1.8);
            plot(x_vec, C_obs_all{sp}(yi, :), '--', ...
                'Color', [0.2 0.2 0.2], 'LineWidth', 1.2, ...
                'Marker', sp_marker{sp}, 'MarkerSize', 3, ...
                'MarkerIndices', 1:max(1, round(Nx/12)):Nx);
            hold off;

            if sp == 1
                title(sprintf('y = %.2f mm', y_vec(yi)), 'FontSize', 9);
            end
            if c == 1
                ylabel(sprintf('%s (%s)', species_names{sp}, unit_C), ...
                       'FontSize', 9, 'Interpreter', 'tex');
            end
            if sp == 3
                xlabel('x (mm)', 'FontSize', 9);
            end
            set(gca, 'FontSize', 8);
            if sp == 1 && c == n_cuts
                legend({'modelled', 'observed'}, 'Location', 'best', 'FontSize', 7);
            end
        end
    end
    sgtitle('Transverse Concentration Profiles  —  modelled vs observed', ...
            'FontSize', 12, 'FontWeight', 'bold');
    save_figure_safe(fig1, 'line_profiles_C_transverse.png', 600);
    fprintf('   Saved -> line_profiles_C_transverse.png\n');

    %% ===================================================================
    %% FIG 2 — Streamwise C(y) at centerline (modelled + observed)
    %% ===================================================================
    fig2 = figure('Name', 'Streamwise Concentration Profiles', ...
                  'Position', [50, 50, 1500, 500], 'Color', 'w');

    subplot(1, 2, 1); hold on; grid on; box on;
    for sp = 1:3
        plot(y_vec, C_th_all{sp}(:, x_idx_center), '-', ...
             'Color', sp_colors{sp}, 'LineWidth', 1.8);
    end
    hold off;
    xlabel('y (mm)'); ylabel(sprintf('C (%s)', unit_C));
    title(sprintf('Modelled  (x = %.3f mm)', x_vec(x_idx_center)));
    legend(species_names, 'Location', 'best', 'Interpreter', 'tex');
    set(gca, 'FontSize', 9);

    subplot(1, 2, 2); hold on; grid on; box on;
    for sp = 1:3
        plot(y_vec, C_obs_all{sp}(:, x_idx_center), '-', ...
             'Color', sp_colors{sp}, 'LineWidth', 1.8);
    end
    hold off;
    xlabel('y (mm)'); ylabel(sprintf('C (%s)', unit_C));
    title(sprintf('Observed  (x = %.3f mm)', x_vec(x_idx_center)));
    legend(species_names, 'Location', 'best', 'Interpreter', 'tex');
    set(gca, 'FontSize', 9);

    sgtitle('Streamwise Concentration along Channel Centerline', ...
            'FontSize', 12, 'FontWeight', 'bold');
    save_figure_safe(fig2, 'line_profiles_C_streamwise.png', 600);
    fprintf('   Saved -> line_profiles_C_streamwise.png\n');

    %% ===================================================================
    %% FIG 3 — Diffusion coefficient line profiles
    %%   Row 1: D(x) transverse at y = 1, 4, 8 mm
    %%   Row 2: D(y) streamwise at centerline (with literature reference)
    %% ===================================================================
    fig3 = figure('Name', 'Diffusion Coefficient Line Profiles', ...
                  'Position', [50, 50, 1500, 850], 'Color', 'w');

    % Transverse sampling stations at fixed physical depths y = 1, 4, 8 mm
    y_station_mm = [1, 4, 8];
    y_stations   = zeros(1, numel(y_station_mm));
    for k = 1:numel(y_station_mm)
        [~, y_stations(k)] = min(abs(y_vec - y_station_mm(k)));
    end

    for sp = 1:3
        subplot(2, 3, sp); hold on; grid on; box on;
        styles = {'-', '--', ':'};
        for k = 1:3
            yi = y_stations(k);
            plot(x_vec, D_all{sp}(yi, :), styles{k}, ...
                 'Color', sp_colors{sp}, 'LineWidth', 1.6);
        end
        yline(D_lit_array(sp), 'k-.', 'LineWidth', 1.2);
        hold off;
        xlabel('x (mm)'); ylabel(sprintf('D (%s)', unit_D));
        title(species_names{sp}, 'Interpreter', 'tex', 'FontSize', 9);
        if sp == 1
            legend({sprintf('y = %.0f mm', y_vec(y_stations(1))), ...
                    sprintf('y = %.0f mm', y_vec(y_stations(2))), ...
                    sprintf('y = %.0f mm', y_vec(y_stations(3))), ...
                    'literature'}, 'Location', 'best', 'FontSize', 7);
        end
        set(gca, 'FontSize', 8);

        subplot(2, 3, 3 + sp); hold on; grid on; box on;
        plot(y_vec, D_all{sp}(:, x_idx_center), '-', ...
             'Color', sp_colors{sp}, 'LineWidth', 1.8);
        yline(D_lit_array(sp), 'k-.', 'LineWidth', 1.2);
        ylim_now = ylim;
        plot([4 4], ylim_now, 'r--', 'LineWidth', 1.2);
        hold off;
        xlabel('y (mm)'); ylabel(sprintf('D (%s)', unit_D));
        title(sprintf('%s  centerline', species_names{sp}), ...
              'Interpreter', 'tex', 'FontSize', 9);
        if sp == 1
            legend({'estimated', 'literature', 'y = 4 mm'}, ...
                   'Location', 'best', 'FontSize', 7);
        end
        set(gca, 'FontSize', 8);
    end

    sgtitle('Diffusion Coefficient Line Profiles  (with literature reference)', ...
            'FontSize', 12, 'FontWeight', 'bold');
    save_figure_safe(fig3, 'line_profiles_D_lines.png', 300);
    fprintf('   Saved -> line_profiles_D_lines.png\n');

    %% ===================================================================
    %% FIG 4 — Source term line profiles
    %%   Row 1: S(x) transverse at y = 1, 4, 8 mm
    %%   Row 2: S(y) streamwise at centerline
    %% ===================================================================
    fig4 = figure('Name', 'Source Term Line Profiles', ...
                  'Position', [50, 50, 1500, 850], 'Color', 'w');

    for sp = 1:3
        subplot(2, 3, sp); hold on; grid on; box on;
        styles = {'-', '--', ':'};
        for k = 1:3
            yi = y_stations(k);
            plot(x_vec, S_all{sp}(yi, :), styles{k}, ...
                 'Color', sp_colors{sp}, 'LineWidth', 1.6);
        end
        yline(0, 'k:', 'LineWidth', 0.8);
        hold off;
        xlabel('x (mm)'); ylabel(sprintf('S (%s)', unit_S));
        title(species_names{sp}, 'Interpreter', 'tex', 'FontSize', 9);
        if sp == 1
            legend({sprintf('y = %.0f mm', y_vec(y_stations(1))), ...
                    sprintf('y = %.0f mm', y_vec(y_stations(2))), ...
                    sprintf('y = %.0f mm', y_vec(y_stations(3)))}, ...
                   'Location', 'best', 'FontSize', 7);
        end
        set(gca, 'FontSize', 8);

        subplot(2, 3, 3 + sp); hold on; grid on; box on;
        plot(y_vec, S_all{sp}(:, x_idx_center), '-', ...
             'Color', sp_colors{sp}, 'LineWidth', 1.8);
        yline(0, 'k:', 'LineWidth', 0.8);
        ylim_now = ylim;
        plot([4 4], ylim_now, 'r--', 'LineWidth', 1.2);
        hold off;
        xlabel('y (mm)'); ylabel(sprintf('S (%s)', unit_S));
        title(sprintf('%s  centerline', species_names{sp}), ...
              'Interpreter', 'tex', 'FontSize', 9);
        set(gca, 'FontSize', 8);
    end

    sgtitle('Source Term Line Profiles  (Gaussian-distributed)', ...
            'FontSize', 12, 'FontWeight', 'bold');
    save_figure_safe(fig4, 'line_profiles_S_lines.png', 600);
    fprintf('   Saved -> line_profiles_S_lines.png\n');

    %% ===================================================================
    %% FIG 5 — Initial estimates D(x) and S(x) at y = 1, 4, 8 mm (Fig. 4 of the article)
    %%   Column (a): D(x); column (b): S(x); one row per species; 600 dpi export.
    %% ===================================================================
    y_section_mm = [1.0, 4.0, 8.0];
    n_st = numel(y_section_mm);
    y_st = zeros(1, n_st);
    for i = 1:n_st
        [~, y_st(i)] = min(abs(y_vec - y_section_mm(i)));
    end
    st_styles = {'-', '--', '-.'};
    xlims = [min(x_vec), max(x_vec)];

    fig5 = figure('Name', 'D-S Cross-Section Profiles (Fig. 4 of the article)', ...
                  'Position', [50, 50, 1150, 900], 'Color', 'w');

    for sp = 1:3
        % --- Column 1 = panel (a): diffusion coefficient D(x) ------------
        subplot(3, 2, (sp-1)*2 + 1); hold on; grid on; box on;
        for k = 1:n_st
            plot(x_vec, D_all{sp}(y_st(k), :), st_styles{k}, ...
                 'Color', sp_colors{sp}, 'LineWidth', 1.6);
        end
        yline(D_lit_array(sp), 'k-.', 'LineWidth', 1.0);
        hold off;
        xlim(xlims); yl = ylim; ylim([0, yl(2)]);
        ylabel(sprintf('D_{%s} (%s)', species_tags{sp}, unit_D), ...
               'Interpreter', 'tex');
        if sp == 3, xlabel('x (mm)'); end
        if sp == 1
            title('', 'FontSize', 11, 'FontWeight', 'bold');
            legend({'y = 1 mm', 'y = 4 mm', 'y = 8 mm', 'literature'}, ...
                   'Location', 'best', 'FontSize', 8, 'Box', 'off');
        end
        set(gca, 'FontSize', 9);

        % --- Column 2 = panel (b): source term S(x) ---------------------
        subplot(3, 2, (sp-1)*2 + 2); hold on; grid on; box on;
        for k = 1:n_st
            plot(x_vec, S_all{sp}(y_st(k), :), st_styles{k}, ...
                 'Color', sp_colors{sp}, 'LineWidth', 1.6);
        end
        yline(0, 'k:', 'LineWidth', 0.8);
        hold off;
        xlim(xlims);
        ylabel(sprintf('S_{%s} (%s)', species_tags{sp}, unit_S), ...
               'Interpreter', 'tex');
        if sp == 3, xlabel('x (mm)'); end
        if sp == 1
            title('', 'FontSize', 11, 'FontWeight', 'bold');
        end
        set(gca, 'FontSize', 9);
    end

    save_figure_safe(fig5, 'line_profiles_DS.png', 600);
    fprintf('   Saved -> line_profiles_DS.png\n');
end

%% ========================================================================
%% LOCAL FUNCTION: 3D SURFACE PLOTS OF D AND S
%% ========================================================================

function create_3d_profiles_DS(X, Y, D_all, S_all, species_names, D_lit_array, valid_masks, opts)
% CREATE_3D_PROFILES_DS  3D surfaces of D(x,y) and S(x,y) for the three species.
%
%   X, Y         : Ny x Nx physical grids (mm)
%   D_all        : 1x3 cell, D in mm^2/s (NaN where not estimated)
%   S_all        : 1x3 cell, S in mol/(L s)
%   species_names: 1x3 cell of names (TeX)
%   D_lit_array  : 1x3 literature D (mm^2/s)
%   valid_masks  : 1x3 cell of logical masks (used only if opts.mask = true)
%   opts         : mask, robust, stations, view, dpi (see Section 7e)
%
%   D is plotted in units of 1e-3 mm^2/s so the z ticks stay readable.
%   Output: profiles_3D_DS.png

    if nargin < 8, opts = struct(); end
    if ~isfield(opts, 'mask'),     opts.mask     = false;       end
    if ~isfield(opts, 'robust'),   opts.robust   = true;        end
    if ~isfield(opts, 'stations'), opts.stations = [1 4 8];     end
    if ~isfield(opts, 'view'),     opts.view     = [-37.5 28];  end
    if ~isfield(opts, 'dpi'),      opts.dpi      = 600;         end

    x  = X(1, :);  y = Y(:, 1);
    xl = [min(x) max(x)];  yl = [min(y) max(y)];
    [XP, YP] = meshgrid(xl, yl);          % corners of the reference plane
    letters = 'abcdef';
    Dscale  = 1e3;                        % D shown as x10^-3 mm^2/s

    st_idx = zeros(1, numel(opts.stations));
    for k = 1:numel(opts.stations)
        [~, st_idx(k)] = min(abs(y - opts.stations(k)));
    end

    fig = figure('Name', '3D Profiles of D and S', ...
                 'Position', [60, 60, 1500, 860], 'Color', 'w');
    tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    for row = 1:2
        for sp = 1:3
            ax = nexttile(tl);
            panel = sprintf('(%s)  %s', letters((row-1)*3 + sp), species_names{sp});

            if row == 1
                Z    = D_all{sp} * Dscale;
                ref  = D_lit_array(sp) * Dscale;
                zlab = 'D (\times10^{-3} mm^2/s)';
                cmap = parula(256);
            else
                Z    = S_all{sp};
                ref  = 0;
                zlab = 'S (mol/(L\cdot s))';
                cmap = diverging_map_(256);
            end
            if opts.mask && ~isempty(valid_masks)
                Z(~valid_masks{sp}) = NaN;
            end

            v = Z(isfinite(Z));
            if isempty(v)
                title(ax, [panel '  (no valid data)'], 'FontSize', 11);
                axis(ax, 'off');
                continue;
            end

            % limits (robust = ignore the extreme 1 % at each end)
            if opts.robust
                lo = pct3d_(v, 1);  hi = pct3d_(v, 99);
            else
                lo = min(v);        hi = max(v);
            end
            if row == 1
                cl = [lo hi];
            else
                m  = max(abs([lo hi]));
                cl = [-m m];                  % symmetric: white = 0
            end
            if cl(2) <= cl(1)
                cl = cl(1) + [-1 1] * max(abs(cl(1)) * 0.05, eps);
            end
            zr  = [min([cl(1), ref]), max([cl(2), ref])];
            pad = 0.05 * diff(zr);

            % surface
            surf(ax, X, Y, Z, 'EdgeColor', 'none', 'FaceColor', 'interp');
            hold(ax, 'on');

            % reference plane (literature D, or S = 0)
            surf(ax, XP, YP, ref * ones(2), ...
                 'FaceColor', [0.35 0.35 0.35], 'FaceAlpha', 0.15, ...
                 'EdgeColor', [0.25 0.25 0.25], 'LineStyle', '--');

            % sampling stations
            for k = 1:numel(st_idx)
                plot3(ax, x, y(st_idx(k)) * ones(size(x)), Z(st_idx(k), :), ...
                      'k-', 'LineWidth', 1.3);
            end
            hold(ax, 'off');

            colormap(ax, cmap);  caxis(ax, cl);
            xlim(ax, xl);  ylim(ax, yl);  zlim(ax, zr + [-pad pad]);
            view(ax, opts.view);  pbaspect(ax, [1 2 0.9]);
            grid(ax, 'on');  box(ax, 'on');
            set(ax, 'FontSize', 9);

            xlabel(ax, 'x (mm)');  ylabel(ax, 'y (mm)');
            zlabel(ax, zlab, 'Interpreter', 'tex');
            title(ax, panel, 'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'tex');

            cb = colorbar(ax);
            cb.Label.String = zlab;  cb.Label.Interpreter = 'tex';  cb.FontSize = 8;
        end
    end

    save_figure_safe(fig, 'profiles_3D_DS.png', opts.dpi);
end

function c = diverging_map_(n)
% Blue - white - red map so that S = 0 is white.
    b = [0.23 0.30 0.75];  w = [0.97 0.97 0.97];  r = [0.71 0.02 0.15];
    h = floor(n/2);
    c = [interp1([0 1], [b; w], linspace(0, 1, h)); ...
         interp1([0 1], [w; r], linspace(0, 1, n - h))];
end

function v = pct3d_(x, p)
% Simple percentile (no toolbox needed).
    x = sort(x(:));
    v = x(max(1, min(numel(x), round(p/100 * numel(x)))));
end

%% ========================================================================
%% LOCAL FUNCTION: LEAST-SQUARES PROFILE FIT
%% ========================================================================

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

%% ========================================================================
%% LOCAL FUNCTION: FIGURE EXPORT
%% ========================================================================

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
