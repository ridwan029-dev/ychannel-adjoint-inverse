function K = reg_matrix(Ny, Nx, dx, dy)
% REG_MATRIX  K = Gx'Gx + Gy'Gy with forward-difference gradients (column-major).
% Then  0.5*sum(|grad f|^2)*dx*dy = 0.5*f(:)'*K*f(:)*dx*dy   (H1 seminorm)
% and its L2 gradient is K*f(:)  (a Neumann Laplacian with a minus sign).
    ex = ones(Nx, 1);  ey = ones(Ny, 1);
    Dx1 = spdiags([-ex ex], [0 1], Nx-1, Nx) / dx;   % (Nx-1) x Nx
    Dy1 = spdiags([-ey ey], [0 1], Ny-1, Ny) / dy;   % (Ny-1) x Ny
    Gx = kron(Dx1, speye(Ny));     % derivative along x (slow index)
    Gy = kron(speye(Nx), Dy1);     % derivative along y (fast index)
    K = Gx.'*Gx + Gy.'*Gy;
end
