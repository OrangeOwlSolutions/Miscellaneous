% --- From non-uniform to uniform

clc
clear all
close all

lambda = 1;
beta = 2 * pi / lambda;

c = 1.5;                                                    % --- Oversampling factor >=1 (typically 2)
K = 1;                                                      % --- 2K+1 interpolation samples (N should be >> K)

Tab = [1.5 1 1; 1.5 2 1.1; 2 3 1.2; 2 6 1.1; 1.5 3 1.5; 1.5 6 1.1];

for k = 1 : 4,
    if ((c == Tab(k, 1)) && (K == Tab(k, 2)))
        SBP_factor = Tab(k, 3);
    end
end

SB_Product                      = SBP_factor * ((2 * pi - pi / c) * K);
Max_Num_PFs                     = ceil(2 * SB_Product / pi);                          % --- Maximum number of PFs

% --- Number of (even) array elements
% N = 20;
N = 14;

% --- Output spectral points
M = 20 * N;
u = 2 * beta * (-M / 2 : (M / 2 - 1)) / M;

% --- Array element locations
% x = [0.139 0.649 0.979 1.501 1.934 2.576 3.261] * lambda;
x = [0.25 0.75 1.25 1.75 2.25 2.75 3.25] * lambda;
x = [-x x];

data = ones(size(x));

result_NFFT_BLAS                = NFFT2_1D_BLAS(data, 2 * beta * x / M, M);
result_NFFT_Matlab              = NFFT2_1D(data, 2 * x / lambda, M, c, K);
result_Gaussian_NFFT_Matlab     = NFFT2_Gaussian_1D(data, 2 * x / lambda, M, c, K);
result_Opt_NUFFT_Matlab         = NFFT2_1D_Opt(data, 2 * x / lambda, M, c, K, Max_Num_PFs, SBP_factor);

rms_NFFT                        = 100*sqrt(sum(abs(result_NFFT_Matlab           - result_NFFT_BLAS.') .^2 ) / sum(abs(result_NFFT_BLAS) .^ 2));
rms_Gaussian_NFFT               = 100*sqrt(sum(abs(result_Gaussian_NFFT_Matlab  - result_NFFT_BLAS.') .^2 ) / sum(abs(result_NFFT_BLAS) .^ 2));
rms_Opt_NFFT                    = 100*sqrt(sum(abs(result_Opt_NUFFT_Matlab      - result_NFFT_BLAS.') .^2 ) / sum(abs(result_NFFT_BLAS) .^ 2));


figure(1)
plot(u/beta,20*log10(abs(result_NFFT_BLAS)),'LineWidth',2)
hold on
plot(u/beta,20*log10(abs(result_NFFT_Matlab)),'r-.','LineWidth',2)
hold off
axis([-1 1 -15 25])
xlabel('u/\beta')
ylabel('20log_{10}(f)')
legend('Exact','Kaiser-Bessel')
title('Kaiser-Bessel')
set(gca,'FontSize',13)
set(findall(gcf,'type','text'),'FontSize',13)

figure(2)
plot(u/beta,20*log10(abs(result_NFFT_BLAS)),'LineWidth',2)
hold on
plot(u/beta,20*log10(abs(result_Gaussian_NFFT_Matlab)),'r-.','LineWidth',2)
hold off
axis([-1 1 -15 25])
xlabel('u/\beta')
ylabel('20log_{10}(f)')
legend('Exact','Gaussian')
title('Gaussian')
set(gca,'FontSize',13)
set(findall(gcf,'type','text'),'FontSize',13)

figure(3)
plot(u/beta,20*log10(abs(result_NFFT_BLAS)),'LineWidth',2)
hold on
plot(u/beta,20*log10(abs(result_Opt_NUFFT_Matlab)),'r-.','LineWidth',2)
hold off
axis([-1 1 -15 25])
xlabel('u/\beta')
ylabel('20log_{10}(f)')
legend('Exact','Optimized')
title('Optimized')
set(gca,'FontSize',13)
set(findall(gcf,'type','text'),'FontSize',13)


