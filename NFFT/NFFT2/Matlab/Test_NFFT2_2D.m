clear all
close all
clc

% --- M_x * M_y input points
M_x = 11;
M_y = 11;

% --- N x M output points
N = 22;
M = 22;

lambda = 1;

x = M_x * (lambda / 2) * (rand(1, M_x * M_y) - 0.5);
y = M_y * (lambda / 2) * (rand(1, M_x * M_y) - 0.5);

data = (rand(1, M_x * M_y) - 0.5) + 1i * (rand(1, M_x * M_y) - 0.5);

Max_Num_PFs = 20;
SBP_factor  = 1.1;
c = 2;
K = 6;
transf_1 = NDFT2_2D(data, x, y, N, M);
transf_2 = NFFT2_2D(data, y, x, M, N);
transf_3 = NFFT2_2D_Opt(data, y, x, M, N, c, K, Max_Num_PFs, SBP_factor);

100 * sqrt(sum(sum(abs(transf_1 - transf_2).^2)) / sum(sum(abs(transf_1).^2)))
100 * sqrt(sum(sum(abs(transf_1 - transf_3).^2)) / sum(sum(abs(transf_1).^2)))

figure(1)
imagesc(abs(transf_1)), colorbar
figure(2)
imagesc(abs(transf_2)), colorbar

figure(3)
imagesc(angle(transf_1)), colorbar
figure(4)
imagesc(angle(transf_2)), colorbar
