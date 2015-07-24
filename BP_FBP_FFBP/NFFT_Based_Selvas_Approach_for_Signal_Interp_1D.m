%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NFFT BASED SELVA'S APPROACH TO BANDLIMITED SIGNAL INTERPOLATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sUp = NFFT_Based_Selvas_Approach_for_Signal_Interp_1D(s, T, a, G, kg, xn)

N = length(s);
sF = fft(s);                                % --- Compute the FFT of the signal samples (see eq. (25)) 
sF = sF(1+mod(-kg:kg,N)) .* G * (N*a);      % --- Repeat some of the spectral samples periodically, and apply the window specified by the vector G.
sF = sF(:).';
 
if mod(a*N,2) == 0                          % --- Pad with zeros to obtain over-sampling factor a.
  sF = [zeros(1,N*a/2-kg),sF,zeros(1,N*a/2-kg-1)];
else
  sF = [zeros(1,(N*a-2*kg-1)/2),sF,zeros(1,(N*a-2*kg-1)/2)];
end

sUp = (1 / (a * N)) * conj(NFFT1_1D(conj(sF), xn / (T / a))); % --- Compute the IDFT by a NER NUFFT
