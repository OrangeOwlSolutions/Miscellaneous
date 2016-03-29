clear all
close all
clc

% load Timings_template_double_1GPU_tol10e_7.txt
% timings1 = Timings_template_double_1GPU_tol10e_7;
% 
% load Timing_Double.txt
% timings2 = Timing_Double;

load Timings_template_float_1GPU_tol10e_4.txt
timings1 = Timings_template_float_1GPU_tol10e_4;

load Timing_Single.txt
timings2 = Timing_Single;

kk = 2.^(4 : 20);
mm = 2 : 8;

figure(1)
hold on
col = hsv(length(mm));
h = zeros(1, length(mm));
% --- column #11:   SVD computation without plan, but with memory transfers
columnIndex = 11;
for p = 0 : length(mm) - 1,
    
   A1 = timings1((1 + p * length(kk)): (1 + p * length(kk)) + 17 - 1, :);
   A2 = timings2((1 + p * length(kk)): (1 + p * length(kk)) + 17 - 1, :);
   temp = A2(:, 5) ./ A1(:, columnIndex);
   h(p + 1) = semilogx(kk, temp, 'color', col(p + 1, :), 'LineWidth', 2, 'DisplayName', sprintf('(m x n) = (%d x %d)', mm(p + 1), mm(p + 1)));
    
end
legend(h)
hold off