function PSI = computePSINn(c, Nmax, nmax, x, y)

dx          = x(2) - x(1);
dy          = y(2) - y(1);
[X, Y]      = meshgrid(x, y);
Support      = zeros(size(X));
Support((X.^2 + Y.^2) <= 1) = 1;

Nx          = length(x);
Ny          = length(y);

[THETA, R]  = cart2pol(X, Y);
r           = reshape(R, 1, numel(R));
theta       = reshape(THETA, 1, numel(THETA));

PSI     = zeros(Nmax, nmax, Nx, Ny);
for N = 0 : Nmax - 1,
    RNn = zeros(nmax, length(r));
    ind = find(r <= 1);
    RNn(:, ind) = computeRNn(nmax, N, c, r(ind));
    if N == 0
        for n = 0 : nmax - 1,
            tempr = RNn(n + 1, :) / sqrt(2 * pi);
            TEMP = Support .* reshape(tempr, Nx, Ny) .* cos(N * THETA);
            PSI(N + 1, n + 1, :, :) = TEMP / sqrt((sum(sum(abs(TEMP).^2)) * (dx * dy)));
        end
    elseif (mod(N, 2) == 0)
        for n = 0 : nmax - 1,
            tempr = RNn(n + 1, :) / sqrt(pi);
            TEMP = Support .* reshape(tempr, Nx, Ny) .* cos(N / 2 * THETA);
            PSI(N + 1, n + 1, :, :) = TEMP / sqrt((sum(sum(abs(TEMP).^2)) * (dx * dy)));
        end
    else
        for n = 0 : nmax - 1,
            tempr = RNn(n + 1, :) / sqrt(pi);
            TEMP = Support .* reshape(tempr, Nx, Ny) .* sin((N + 1) / 2 * THETA);
            PSI(N + 1, n + 1, :, :) = TEMP / sqrt((sum(sum(abs(TEMP).^2)) * (dx * dy)));
        end
    end
end

