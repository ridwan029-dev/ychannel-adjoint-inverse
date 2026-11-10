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
