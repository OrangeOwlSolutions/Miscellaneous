clear all
close all
clc

load Timings_template_double_1GPU_tol10e_7.txt
timings = Timings_template_double_1GPU_tol10e_7;

% load Timings_template_double_1GPU_tol_10e_7_m_different_from_n.txt
% timings = Timings_template_double_1GPU_tol_10e_7_m_different_from_n;

kk = 2.^(4 : 20);
mm = 2 : 8;

figure(1)
hold on
col = hsv(length(mm));
h = zeros(1, length(mm));
% --- column #11:   SVD computation without plan, but with memory transfers
columnIndex = 11;
for p = 0 : length(mm) - 1,
    
   A = timings((1 + p * length(kk)): (1 + p * length(kk)) + 17 - 1, :);
   h(p + 1) = semilogx(kk, A(:, columnIndex), 'color', col(p + 1, :), 'LineWidth', 2, 'DisplayName', sprintf('(m x n) = (%d x %d)', mm(p + 1), mm(p + 1)));
    
end
legend(h)
hold off