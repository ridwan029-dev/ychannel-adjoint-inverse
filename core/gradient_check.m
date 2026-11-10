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
