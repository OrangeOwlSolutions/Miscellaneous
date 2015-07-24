%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FREQUENCY DOMAIN INTERPOLATION WITH RECTANGULAR WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = rectNFFT1D(x, y, xn, B)

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
len=length(xn);
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

% --- Differential frequency
df = 1/(N*T);

% --- kg factor
kg = floor((1/T-B/2)/df);  

%--- NUFFT based Selva's approach to bandlimited signal interpolation with
%rectangular window
% result = NUFFT_Based_Selvas_Approach_for_Signal_Interpolation_Rect_1D(y, T, a, kg, xn);
result = NFFT_Based_Selvas_Approach_for_Signal_Interp_1D(y, T, a, 1, kg, xn);
result=result(1:len);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUFFT BASED SELVA'S APPROACH TO BANDLIMITED SIGNAL INTERPOLATION WITH RECTANGULAR WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function sUp = NUFFT_Based_Selvas_Approach_for_Signal_Interpolation_Rect_1D(s, T, a, kg, xn)
% 
% N = length(s);
% sF = fft(s);                                % --- Compute the FFT of the signal samples (see eq. (25)) 
% sF = sF(1+mod(-kg:kg,N))*a;                 % Repeat some of the spectral samples periodically, and apply the window specified by the vector G.
% sF = sF(:).';
%  
% if mod(a*N,2) == 0                          % --- Pad with zeros to obtain over-sampling factor a.
%   sF = [zeros(1,N*a/2-kg),sF,zeros(1,N*a/2-kg-1)];
% else
%   sF = [zeros(1,(N*a-2*kg-1)/2),sF,zeros(1,(N*a-2*kg-1)/2)];
% end
% 
% sUp = (1/(a*N))*conj(NufftNer1D(conj(sF),xn/(T/a))); % --- Compute the IDFT by a NER NUFFT

