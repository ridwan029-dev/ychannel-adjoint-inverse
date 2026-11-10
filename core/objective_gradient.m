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
