function z = det_randn(sz, seed)
%DET_RANDN  Standard normal numbers from a fixed seed (Park-Miller generator
%   and Box-Muller transform), so that the synthetic noise of the tests does
%   not depend on the MATLAB version or its random number settings.
%   z = det_randn([m n], seed) returns an m x n array.
    n = prod(sz);  m = 2147483647;  s = mod(seed*104729 + 12345, m - 1) + 1;
    for k = 1:50, s = mod(16807*s, m); end      % discard the first draws
    u = zeros(2*ceil(n/2), 1);
    for k = 1:numel(u)
        s = mod(16807*s, m);
        u(k) = s/m;
    end
    u1 = u(1:2:end);  u2 = u(2:2:end);
    r = sqrt(-2*log(u1));
    z = [r.*cos(2*pi*u2); r.*sin(2*pi*u2)];
    z = reshape(z(1:n), sz);
end
