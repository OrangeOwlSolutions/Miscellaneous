%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME DOMAIN INTERPOLATION WITH APPROXIMATE PROLATE SPHEROIDAL WAVEFUNCTION WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = knab1D(x, y, xn, B)

% --- Mean sampling step
T = mean(diff(x));
if T > 1 / B     
    error('Nyquist frequency is bigger than sampling frequency');    
end

P = 18;                                 % --- 2*P+1 is the number of samples involved in the sampling series
Tw = P * T;

[TN, X] = meshgrid(xn, x);

% --- Time domain approximate prolate spheroidal wave function
wap = sinc((1 / T - B) * sqrt((TN - X).^2 - Tw^2)) / sinc(1i * (1 / T - B) * Tw);

% --- Time domain ectangular window
g0=sinc(1/T*(TN-X));

result = y * (wap .* g0);
