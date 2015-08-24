function result = NFFT1_1D_BLAS(data, x)

% --- Calculates

% \sum_{n} d_n * exp(-j * k * x_n)

% where x_n are random points, d_n are the data and k = -N / 2 : N / 2 - 1

N=length(data);

k=-N/2:N/2-1;

[K,X]=meshgrid(k,x);

result = exp(-j*K.*X)*data.';
