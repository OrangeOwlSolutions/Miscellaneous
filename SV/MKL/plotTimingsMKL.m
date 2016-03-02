clear all
close all
clc

load Timing_Double.txt
timings = Timing_Double;

kk = 2.^(4 : 20);
mm = 2 : 8;

figure(1)
hold on
col = hsv(length(mm));
h = zeros(1, length(mm));
for p = 0 : length(mm) - 1,
    
   A = timings((1 + p * length(kk)): (1 + p * length(kk)) + 17 - 1, :);
   h(p + 1) = semilogx(kk, A(:, 5), 'color', col(p + 1, :), 'LineWidth', 2, 'DisplayName', sprintf('(m x n) = (%d x %d)', mm(p + 1), mm(p + 1)));
    
end
legend(h)
hold off
