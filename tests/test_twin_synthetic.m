%TEST_TWIN_SYNTHETIC  Synthetic twin test (Section 6.5): known D and S ->
%   synthetic data -> inversion, on the grid of the experiment
%   (270 x 81, x = -0.5..0.5 mm, y = 0..8.07 mm, v = 0.15 mm/s).
%   Checks: (1) exact gradient, (2) noise-free data, four CG methods: decrease
%   of J and errors of D and S (identifiability), (3) noisy data (0.03 mol/L): the four CG
%   methods and the iteration at which the data misfit reaches the noise level.
%   Output: twin_results.mat in this folder. Uses ../core; base MATLAB only.
%   Part of ychannel-adjoint-inverse; see README.md for the article reference.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'core'));  addpath(here);

Nx = 81; Ny = 270;
x = linspace(-0.5, 0.5, Nx);  y = (0:Ny-1).'*0.03;     % mm, grid of the experiment
[X, Y] = meshgrid(x, y);
v = 0.15;                                     % mm/s

% ---- true fields ----
Db = [1.76e-3, 2.03e-3, 0.81e-3];            % mm^2/s, reference values at 28 degC (Table 2)
Gd = exp(-X.^2/(2*0.03^2));
Gs = exp(-X.^2/(2*0.025^2));
a  = 0.05*exp(-Y/4);                          % mol/(L s)
Dt = {Db(1)*(1 + 1.5*Gd), Db(2)*(1 + 1.5*Gd), Db(3)*(1 + 1.5*Gd)};
St = {-a.*Gs, -2*a.*Gs, a.*Gs};

% ---- boundary values: inlet profiles, walls = inlet wall values ----
Cin = {1.0*0.5*erfc(X/0.01), 2.0*0.5*erfc(-X/0.01), zeros(Ny, Nx)};
opt = struct('alphaD', 1e-7, 'alphaS', 1e-7, 'Dref', 1e-3, 'Sref', 1e-2, ...
             'Dmin', 1e-5, 'Dmax', 1e-1);
W = {ones(Ny, Nx), ones(Ny, Nx), ones(Ny, Nx)};
prob0 = make_problem(x, y, v, Cin, W, Cin, opt);
Ct = cell(1, 3);
for s = 1:3, Ct{s} = forward_solve(Dt{s}, St{s}, s, prob0); end

% ---- initial guess: bulk D everywhere, S from the integral balance ----
D0 = {Db(1)*ones(Ny, Nx), Db(2)*ones(Ny, Nx), Db(3)*ones(Ny, Nx)};
S0 = cell(1, 3);
for s = 1:3, S0{s} = init_source_balance(Ct{s}, x, y, v, D0{s}, 0.03); end
band = abs(X) < 0.08 & Y > 0.5;
relerr = @(A, B) norm(A(band) - B(band))/norm(B(band));
errD = @(r) mean([relerr(r.D1, Dt{1}), relerr(r.D2, Dt{2}), relerr(r.D3, Dt{3})]);
errS = @(r) mean([relerr(r.S1, St{1}), relerr(r.S2, St{2}), relerr(r.S3, St{3})]);

% ---- (1) gradient check ----
prob = make_problem(x, y, v, Ct, W, Ct, opt);
fprintf('\n--- gradient check at the initial guess ---\n');
gradient_check(D0, S0, prob, 2);

% ---- (2) noise-free data, all four CG methods ----
methods = {'FR', 'PR', 'PB', 'DY'};
fprintf('\n--- noise-free data (300 iterations) ---\n');
fprintf('%-4s %-18s %6s %11s %11s %7s %7s\n', 'CG', 'stop', 'iter', 'J0', 'J_end', 'errD', 'errS');
res0 = cell(1, numel(methods));
for m = 1:numel(methods)
    r = cg_optimize(prob, D0, S0, methods{m}, struct('maxIter', 300, 'verbose', false));
    res0{m} = r;
    fprintf('%-4s %-18s %6d %11.4e %11.4e %7.3f %7.3f\n', methods{m}, r.stop_reason, r.num_updates, ...
        r.initial_cost, r.final_cost, errD(r), errS(r));
end

% ---- (3) noisy data ----
noise = 0.03;                                 % mol/L
Cobs = cell(1, 3);
for s = 1:3, Cobs{s} = Ct{s} + noise*det_randn([Ny Nx], 7 + s); end
prob = make_problem(x, y, v, Cobs, W, Ct, opt);
noise_floor = 0.5*3*Nx*Ny*noise^2*prob.dx*prob.dy;
cgopt = struct('maxIter', 100, 'verbose', false, 'Cpde', {Ct});
fprintf('\n--- twin inversion (noise %.2f mol/L) ---\n', noise);
fprintf('expected data misfit at the noise level: %.4e\n', noise_floor);
fprintf('%-4s %-18s %6s %11s %11s %11s %8s %7s %7s\n', 'CG', 'stop', 'iter', 'J0', 'J_end', 'Jdata_end', 'k_noise', 'errD', 'errS');
res = cell(1, numel(methods));
for m = 1:numel(methods)
    r = cg_optimize(prob, D0, S0, methods{m}, cgopt);
    res{m} = r;
    k = find(r.J_data_history <= noise_floor, 1) - 1;    % update at which the noise level is reached
    if isempty(k), k = NaN; end
    fprintf('%-4s %-18s %6d %11.4e %11.4e %11.4e %8g %7.3f %7.3f\n', methods{m}, r.stop_reason, ...
        r.num_updates, r.initial_cost, r.final_cost, r.J_data_history(end), k, errD(r), errS(r));
end
fprintf('initial guess: relative error D %.3f, S %.3f\n', ...
    mean([relerr(D0{1}, Dt{1}), relerr(D0{2}, Dt{2}), relerr(D0{3}, Dt{3})]), ...
    mean([relerr(S0{1}, St{1}), relerr(S0{2}, St{2}), relerr(S0{3}, St{3})]));
save(fullfile(here, 'twin_results.mat'), 'res', 'res0', 'Dt', 'St', 'D0', 'S0', 'X', 'Y', 'x', 'y');
