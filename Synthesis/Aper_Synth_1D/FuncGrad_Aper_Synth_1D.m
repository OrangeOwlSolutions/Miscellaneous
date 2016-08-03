function Func = FuncGrad_Aper_Synth_1D(c)
    
global Ea Nu uMax delta_uu xsi dmin delta_xsi lambda xsi_nodes

% --- Chebyshev
% f0 = InterpChebyshev(c);                                        % --- Update the polynomial coefficients of f0(xsi)
% f = expansiveFuncGen(f0);
% x = dmin / delta_xsi * polyval(f, xsi);

% --- Lagrange
f0  = mapLagrangeInterpFunction(xsi_nodes, c);              % Coefficienti polinomiali della f0(xsi)
f   = expansiveFuncGen(f0);
x   = polyval(f, xsi);

Spectrum = NFFT2_1D(Ea, x / (lambda / 2), Nu, 2, 6);
Projected_Spectrum = Spectral_projection(Spectrum);
 
Func = sum((abs(Spectrum).^2 - abs(Projected_Spectrum).^2).^2) * delta_uu;

end

