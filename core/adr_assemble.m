function A = adr_assemble(D, prob)
% ADR_ASSEMBLE  Sparse matrix of the discrete Eq. (2) for one species:
%
%     v(x) dC/dy = d/dx( D dC/dx ) + S          (steady, depth-averaged)
%
% Discretisation (conservative, finite volume in x, implicit march in y):
%   v_j (C_ij - C_{i-1,j})/dy
%     - [ D_{i,j+1/2}(C_{i,j+1}-C_ij) - D_{i,j-1/2}(C_ij-C_{i,j-1}) ]/dx^2 = S_ij
%   with D_{j+1/2} = (D_j + D_{j+1})/2.
% Rows for boundary nodes (inlet y = 0 and the two walls x = x_min, x_max)
% are identity rows: C = C_bc (Dirichlet values taken from the data).
% Outlet (y = y_max) needs no condition because the problem is parabolic in y.
%
% Indexing: MATLAB column-major, k = i + (j-1)*Ny, i = y-index, j = x-index.

    Ny = prob.Ny;  Nx = prob.Nx;  dx = prob.dx;  dy = prob.dy;
    N  = Ny*Nx;
    idx = reshape(1:N, Ny, Nx);

    v = prob.v;                        % scalar or 1 x Nx
    if isscalar(v), v = v*ones(1, Nx); end
    V = repmat(v(:).', Ny, 1);

    Dp = 0.5*(D(:, 2:Nx-1) + D(:, 3:Nx));     % D_{j+1/2}, interior columns
    Dm = 0.5*(D(:, 2:Nx-1) + D(:, 1:Nx-2));   % D_{j-1/2}
    Dp = Dp(2:Ny, :);  Dm = Dm(2:Ny, :);
    Vi = V(2:Ny, 2:Nx-1);

    kc = idx(2:Ny,   2:Nx-1);   % centre
    ks = idx(1:Ny-1, 2:Nx-1);   % upstream (i-1)
    ke = idx(2:Ny,   3:Nx);     % j+1
    kw = idx(2:Ny,   1:Nx-2);   % j-1

    rows = [kc(:); kc(:); kc(:); kc(:)];
    cols = [kc(:); ks(:); ke(:); kw(:)];
    vals = [Vi(:)/dy + (Dp(:) + Dm(:))/dx^2;
            -Vi(:)/dy;
            -Dp(:)/dx^2;
            -Dm(:)/dx^2];

    kb = idx(prob.bnd);                       % boundary nodes -> identity rows
    rows = [rows; kb(:)];  cols = [cols; kb(:)];  vals = [vals; ones(numel(kb), 1)];

    A = sparse(rows, cols, vals, N, N);
end
