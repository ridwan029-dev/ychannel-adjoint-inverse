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
