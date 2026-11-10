function [C, A] = forward_solve(D, S, s, prob)
% FORWARD_SOLVE  Solve the discrete Eq. (2) for species s.
    A = adr_assemble(D, prob);
    b = S(:);
    b(prob.bnd(:)) = prob.gbc{s};
    C = reshape(A \ b, prob.Ny, prob.Nx);
end
