%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUFFT BASED SELVA'S APPROACH TO BANDLIMITED SIGNAL INTERPOLATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sUp = NFFT_Based_Selvas_Approach_for_Signal_Interp_2D(INPUT_DATA, T1, T2, a, b, G, kgn, kgm, RHO_OUT_SHIFTED, COSALPHA_OUT_SHIFTED)

N = size(INPUT_DATA, 1);
M = size(INPUT_DATA, 2);

sF = fft2(INPUT_DATA);                              % --- Compute the FFT of the signal samples (see eq. (25)) 
sF = sF(1+mod(-kgn:kgn,N),1+mod(-kgm:kgm,M));
sF = sF .* G * (N*a) * (M*b);                       % --- Repeat some of the spectral samples periodically, and apply the window specified by the vector G.
        
csF = size(sF, 1);
if mod(b*M,2) == 0                                  % --- Pad with zeros to obtain over-sampling factor b.
    sF = [zeros(csF,M*b/2-kgm),sF,zeros(csF,M*b/2-kgm-1)];
else
    sF = [zeros(csF,(M*b-2*kgm-1)/2),sF,zeros(csF,(M*b-2*kgm-1)/2)];
end

rsF=size(sF,2);
if mod(a*N,2) == 0                                  % --- Pad with zeros to obtain over-sampling factor a.
    sF = [zeros(N*a/2-kgn,rsF);sF;zeros(N*a/2-kgn-1,rsF)];
else
    sF = [zeros((N*a-2*kgn-1)/2,rsF);sF;zeros((N*a-2*kgn-1)/2,rsF)];
end

% --- Compute the IDFT by a NER NUFFT
sUp = (1/((a*b*N*M)))*conj(NFFT1_2D(conj(sF),reshape(RHO_OUT_SHIFTED/(T1/((a))),1,(numel(RHO_OUT_SHIFTED))),reshape(COSALPHA_OUT_SHIFTED/(T2/((b))),1,(numel(COSALPHA_OUT_SHIFTED))))); 
