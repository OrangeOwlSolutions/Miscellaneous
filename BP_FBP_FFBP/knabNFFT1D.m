%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FREQUENCY DOMAIN INTERPOLATION WITH APPROXIMATE PROLATE SPHEROIDAL WAVEFUNCTION WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = knabNFFT1D(x, y, xn, B)

% --- B:        Bandwidth of the signal to be interpolated
% --- x:        input sampling points
% --- xn:       output sampling points
% --- y:        input samples

% --- Mean sampling step
T = mean(diff(x));
if T > 1 / B     
    error('Nyquist frequency is bigger than sampling frequency');    
end

% --- Shifting the output sampling locations
xn = (xn - min(x));

% --- Determining the number of input samples
len = length(xn);
N = ((max(xn) - min(xn)) / (T * (1 - 1 / len)));
if (N < 1) || (isnan(N))
    N = 1;
end

% --- Oversampling detection
a = round(len / N);
if a < 2
    a = 2;
    ii = round(a * N) - len;
    xn = [xn, zeros(1, ii)];
    
end

% --- Selva's window
P = 18;                                 % --- 2*P+1 is the number of samples involved in the sampling series
[G, kg] = KnabWindow(T, B, P, length(y));

%--- NUFFT based Selva's approach to bandlimited signal interpolation with
%approximate prolate spheroidal wavefunction window
result = NFFT_Based_Selvas_Approach_for_Signal_Interp_1D(y, T, a, G, kg, xn);

result=result(1: len);


