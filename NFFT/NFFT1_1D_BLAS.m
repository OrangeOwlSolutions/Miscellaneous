function result = NFFT1_1D_BLAS(data, x)

% --- Calculates

% \sum_{k} d_k * exp(-j * k * x_l)

% where x_l are random output sampling locations, d_k are the data and k = -N / 2 : N /
% 2 - 1 and the input sampling locations

N=length(data);

k=-N/2:N/2-1;

[K,X]=meshgrid(k,x);

result = exp(-1i*K.*X)*data.';
