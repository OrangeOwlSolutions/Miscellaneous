function result = NFFT2_1D_BLAS(data, x, N)

% --- Calculates

% \sum_{l} d_l * exp(-j * k * x_l)

% where x_l are random input sampling locations, d_l are the data and k = -N / 2 : N /
% 2 - 1 and the output sampling locations

k=-N/2:N/2-1;

[X,K]=meshgrid(x,k);

result = exp(-1i*K.*X)*data.';
