function transf = NDFT2_2D(data, x, y, N, M)

u = -(N / 2) : (N / 2) - 1;
v = -(M / 2) : (M / 2) - 1;

[U, V] = meshgrid(u, v);

U = reshape(U, N * M, 1);
V = reshape(V, N * M, 1);

[X, U] = meshgrid(x, U);
[Y, V] = meshgrid(y, V);

Kernel = exp(-1i * 2 * pi * X .* U / N) .* exp(-1i * 2 * pi * Y .* V / N);

transf = reshape(Kernel * data.', N, M);
