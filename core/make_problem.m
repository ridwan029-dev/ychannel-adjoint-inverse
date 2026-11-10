function prob = make_problem(x, y, v, Cobs, W, Cbc, opt)
% MAKE_PROBLEM  Collect grid, data and settings for the inverse problem.
%   x (1 x Nx, mm), y (Ny x 1 or 1 x Ny, mm), v (mm/s, scalar or 1 x Nx)
%   Cobs{s}  observed concentration used in the data misfit (mol/L), Ny x Nx
%   W{s}     misfit weights (1 = use point, 0 = ignore, e.g. NaN in data)
%   Cbc{s}   field whose inlet row and wall columns give the Dirichlet values
%   opt      struct: alphaD, alphaS, Dref, Sref, Dmin, Dmax

    prob.x = x(:).';  prob.y = y(:);
    prob.Nx = numel(x);  prob.Ny = numel(y);
    prob.dx = x(2) - x(1);  prob.dy = y(2) - y(1);
    prob.v  = v;
    prob.nsp = numel(Cobs);

    bnd = false(prob.Ny, prob.Nx);
    bnd(1, :) = true;          % inlet
    bnd(:, 1) = true;          % wall x = x_min
    bnd(:, end) = true;        % wall x = x_max
    prob.bnd = bnd;
    prob.interior = ~bnd;

    for s = 1:prob.nsp
        c = Cobs{s};
        w = double(W{s});
        w(isnan(c)) = 0;
        c(isnan(c)) = 0;
        prob.Cobs{s} = c;
        prob.W{s} = w;
        g = Cbc{s};
        g(isnan(g)) = 0;
        prob.gbc{s} = g(bnd);   % Dirichlet values on boundary nodes
    end

    prob.alphaD = opt.alphaD;  prob.alphaS = opt.alphaS;
    prob.Dref = opt.Dref;      prob.Sref = opt.Sref;
    prob.Dmin = opt.Dmin;      prob.Dmax = opt.Dmax;
    prob.K = reg_matrix(prob.Ny, prob.Nx, prob.dx, prob.dy);
end
