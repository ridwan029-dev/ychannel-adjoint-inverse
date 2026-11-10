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
