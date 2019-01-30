% --- From non-uniform to uniform

clc
clear all
close all

lambda = 1;
beta = 2 * pi / lambda;

c = 2;                                                    % --- Oversampling factor >=1 (typically 2)
K = 6;                                                      % --- 2K+1 interpolation samples (N should be >> K)

Tab = [2 3 1.2; 2 6 1.1; 1.5 3 1.5; 1.5 6 1.1];

for k = 1 : 4,
    if ((c == Tab(k, 1)) && (K == Tab(k, 2)))
        SBP_factor = Tab(k, 3);
    end
end

SB_Product                      = SBP_factor * ((2 * pi - pi / c) * K);
Max_Num_PFs                     = ceil(2 * SB_Product / pi);                          % --- Maximum number of PFs

% --- Half-size
a = 20 * lambda;

% --- Number of (even) array elements
N = ceil(2 * a / (lambda / 2));
if (mod(N, 2) ~= 0) N = N + 1; end

% --- Output spectral points
M = N;
u = -M / 2 : (M / 2 - 1);

Num_tests = 100;

rms_Opt_NFFT        = zeros(1, Num_tests);
rms_NFFT            = zeros(1, Num_tests);
rms_Gaussian_NFFT   = zeros(1, Num_tests);

err_Opt_NFFT        = zeros(1, Num_tests);
err_NFFT            = zeros(1, Num_tests);
err_Gaussian_NFFT   = zeros(1, Num_tests);

for tt = 1 : Num_tests,

    % --- Array element locations
    x = 2 * a * (rand(1, N) - 0.5);

    data = randn(1, N) + 1i * randn(1, N);

    result_NFFT_BLAS                = NFFT2_1D_BLAS(data, 2 * beta * x / M, M);
    result_NFFT_Matlab              = NFFT2_1D(data, 2 * x / lambda, M, c, K);
    result_Gaussian_NFFT_Matlab     = NFFT2_Gaussian_1D(data, 2 * x / lambda, M, c, K);
    result_Opt_NFFT_Matlab          = NFFT2_1D_Opt(data, 2 * x / lambda, M, c, K, Max_Num_PFs, SBP_factor);

    rms_NFFT(tt)                    = 100*sqrt(sum(abs(result_NFFT_Matlab           - result_NFFT_BLAS.') .^2 ) / sum(abs(result_NFFT_BLAS) .^ 2));
    rms_Gaussian_NFFT(tt)           = 100*sqrt(sum(abs(result_Gaussian_NFFT_Matlab  - result_NFFT_BLAS.') .^2 ) / sum(abs(result_NFFT_BLAS) .^ 2));
    rms_Opt_NFFT(tt)                = 100*sqrt(sum(abs(result_Opt_NFFT_Matlab - result_NFFT_BLAS.') .^ 2) / sum(abs(result_NFFT_BLAS) .^ 2));
    
    err_NFFT(tt)                    = max(abs(result_NFFT_Matlab           - result_NFFT_BLAS.'));
    err_Gaussian_NFFT(tt)           = max(abs(result_Gaussian_NFFT_Matlab  - result_NFFT_BLAS.'));
    err_Opt_NFFT(tt)                = max(abs(result_Opt_NFFT_Matlab       - result_NFFT_BLAS.'));

    tt / Num_tests

end

mean(rms_NFFT)
mean(rms_Gaussian_NFFT)
mean(rms_Opt_NFFT)

max(err_NFFT)
max(err_Gaussian_NFFT)
max(err_Opt_NFFT)

figure(1)
semilogy(rms_NFFT,'LineWidth',2)
hold on
semilogy(rms_Opt_NFFT,'r','LineWidth',2)
semilogy(rms_Gaussian_NFFT,'g','LineWidth',2)
hold off
xlabel('Realization')
ylabel('Percentage rms error')
legend('Kaiser-Bessel','Optimized','Gaussian')

% —- For drawings only
set(gca,'FontSize',13)
set(findall(gcf,'type','text'),'FontSize',13)

