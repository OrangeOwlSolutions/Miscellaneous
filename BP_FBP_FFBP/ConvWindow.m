%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPROXIMATE PROLATE SPHEROIDAL WAVEFUNCTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- See page 5470, J. Selva, Convolution-Based Trigonometric Interpolation of Band-Limited Signals
function W = ConvWindow(f, B, T)

W = zeros(size(f));

I = abs(f) < B/2;

W(I) = besseli(0,pi*B*(T/2)*sqrt(1-(2*f(I)/B).^2))/(B*sinc(1i*B*T/2));

