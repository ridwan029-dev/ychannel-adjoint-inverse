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
