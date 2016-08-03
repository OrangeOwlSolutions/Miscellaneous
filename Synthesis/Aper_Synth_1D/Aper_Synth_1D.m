clc
clear all
close all

global lambda alfa dmin delta_xsi xsi uMax uu delta_uu N Nu C_ord Ea N_nodes xsi_nodes Upper_mask Lower_mask

lambda = 1;                                                     % --- Wavelength

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS OF THE APERIODIC ARRAY %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 25;                                                         % --- Number of elements
dmin = lambda / 3;                                              % --- Minimum allowed interelement spacing
% dmax = 3 * lambda;                                              % --- Maximum allowed interelement spacing
% alfa = dmax / dmin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCRETIZATION OF THE CONTINUOUS APERTURE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aap = 7 * lambda;                                               % --- Half-size of the continuous aperture
Ppr = ceil(aap / (lambda / 20));                                % --- 2 * Ppr + 1 is the number of discretization points of the refined aperture 
xap_refined = linspace(-aap, aap, 2 * Ppr + 1);                 % --- Refined aperture sampling
xap = -aap : lambda / 2 : aap;                                  % --- Rough aperture sampling

%%%%%%%%%%%%%%%%%%%%
% MAPPING FUNCTION %
%%%%%%%%%%%%%%%%%%%%
% --- Chebyshev
% xsi = linspace(-1, 1, N);                                       % --- Independent variable of the mapping function
% delta_xsi = 2 / (N - 1);                                        % --- Sampling step in the xsi variable
% C_ord     = 10;                                                 % --- Order of teh Chebyshev interpolation polynomials
% N_nodes   = 10;                                                 % --- Number of Chebyshev interpolation nodes (MUST BE > 3)
% xsi_nodes = cos(pi * ((1 : N_nodes) - 0.5) / N_nodes);            % --- Chebyshev interpolation nodes in the xsi variable

% --- Lagrange
xsi = (0 : N - 1) * dmin;                                         % --- Independent variable of the mapping function
L_ord = 12;                                                       % --- Number of Lagrange interpolation nodes (MUST BE > 3)
xsi_nodes = linspace(0, max(xsi), L_ord);                         % --- Lagrange interpolation nodes in the xsi variable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAGRANGE INTERPOLATION POLYNOMIALS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi     = zeros(L_ord, L_ord);
phider1 = zeros(L_ord, L_ord);
phider2 = zeros(L_ord, L_ord);
phider3 = zeros(L_ord, L_ord);
for k = 1 : L_ord,
    phi(k, :)               = poly(xsi_nodes(2 : L_ord)) ./ prod((xsi_nodes(1) - xsi_nodes(2 : L_ord)));
    phider1(k, 2 : L_ord)   = polyder(phi(k, :));
    phider2(k, 3 : L_ord)   = polyder(phider1(k, :));
    phider3(k, 4 : L_ord)   = polyder(phider2(k, :));
    xsi_nodes               = circshift(xsi_nodes.', -1).';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRUM DISCRETIZATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
chi = 8;                                                        % --- Spectrum oversampling factor
K = 1;                                                          % --- Spectral coverage factor (K = 1 means the whole visible)
uMax = K * pi;                                                  % --- Maximum value of the spectral variable u considered by the synthesis algorithm
Nu = round(chi * 2 * uMax / ((2 * pi / (2 * aap)) * lambda / 2));     
                                                                % --- Number of spectral samples
uu = K * ((2 * pi * (0 : Nu - 1) * (1 / Nu)) - pi);             % --- Discretization points of the u axis
delta_uu = uu(2) - uu(1);                                       % --- Spectral sampling step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCE APERTURE FIELD AND ASSOCIATED SPECTRUM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
decibel = 15;                                                   % --- Aperture field tapering
gamma = decibel / 20 * 1 /(aap^2) * log(10);              
Ea = exp(-gamma * xap.^2);                                      % --- Tapered aperture field with Gaussian tapering
Ea_refined = exp(-gamma*(xap_refined).^2);                      % --- Tapered aperture field with Gaussian tapering (refined case)

Spectrum = NFFT2_1D(Ea, xap / (lambda / 2), Nu, 2, 6);             % --- Spectrum of the aperture field (x must be normalized to lambda / 2)
Spectrum = Spectrum / max(abs(Spectrum));

%%%%%%%%%%%%%%%%%%%%%%%%%
% UPPER AND LOWER MASKS %
%%%%%%%%%%%%%%%%%%%%%%%%%
u_first  = pi / 10;
u_second = pi / 5;

Internal_coverage       = 1.5;
External_coverage       = -25;
Intermediate_coverage   = -35;

Low_spectrum_level      = -30;

% u_first  = pi / 8;
% u_second = pi / 3;
% 
% Internal_coverage       = 1.5;
% External_coverage       = -25;
% Intermediate_coverage   = -35;
% 
% Low_spectrum_level      = -30;

indices_first_interval = find(uu <= u_first  & uu >= -u_first);
indices_third_interval = find(uu >= u_second | uu <= -u_second);

Upper_mask                         = ones(size(uu)) * max(abs(Spectrum)) * 10.^(Intermediate_coverage / 20);             
Upper_mask(indices_first_interval) = max(abs(Spectrum)) * 10^(Internal_coverage / 20);
Upper_mask(indices_third_interval) = max(abs(Spectrum)) * 10^(External_coverage / 20);
         
indices_low_Spectrum = find((20 * log10(abs(Spectrum)) < Low_spectrum_level));

Lower_mask =10 .^ ((20 * log10(Spectrum) - 3) / 20);
Lower_mask(indices_low_Spectrum) = 0;

%%%%%%%%%%%%%%%%%%
% STARTING POINT %
%%%%%%%%%%%%%%%%%%
Density = max(abs(Ea_refined)) ./ abs(Ea_refined);
% --- Chebyshev
% Density(find(Density > alfa)) = alfa;
% --- Lagrange
IntDensity = sum(Density) * (xap_refined(3) - xap_refined(2));

delta_x_prime = 2 * aap / (N - 1);

% --- Starting positions as a function of the element density calculated from the tapered aperture
x(1) = 0;
for kk = 1 : N - 1    
    temp = find((xap_refined + aap) <= (kk * delta_x_prime));
    % --- Chebyshev
%     x(kk+1) = dmin / (delta_xsi * aap) * sum(Density(temp)) * (xap_refined(3) - xap_refined(2));
    % --- Lagrange
    x(kk+1) = 2 * aap * sum(Density(temp)) * (xap_refined(3) - xap_refined(2)) ./ IntDensity;
end
x_tapered_aperture = x;                                         % --- Positions determined from the tapered aperture

Spectrum_tapered_aperture = NFFT2_1D(ones(size(x)), x / (lambda / 2), Nu, 2, 6); % --- Spectrum tapered aperture
% Spectrum_tapered_aperture = Spectrum_tapered_aperture / max(abs(Spectrum_tapered_aperture));   

% --- Lagrange
c = interp1(xsi, x, xsi_nodes);
f0 = mapLagrangeInterpFunction(xsi_nodes, c);                   % --- Coefficients of f0(xsi)
x = polyval(f0, xsi);

% --- Chebyshev
% y = xap_refined / aap;
% c = interp1(y, Density, xsi_nodes);
% f0 = InterpChebyshev(c);                                        % --- Coefficients of f0(xsi)
% f = polyint(f0);
% x = dmin / delta_xsi * polyval(f, xsi);
% x = x - min(x);

xstart = x;                                                     % --- Starting positions

Ea = ones(1, N);                                                % --- Excitation of the aperiodic array
Spectrum_initial = NFFT2_1D(Ea, x / (lambda / 2), Nu, 2, 6);
rifer = max(abs(Spectrum_initial));
Ea = Ea / rifer;
Spectrum_initial = Spectrum_initial / rifer;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIMIZATION OF THE POSITIONS ONLY %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_old = c;
% load Intermediate_result.mat
[cnew, Functional_value, exitFlag, output] = fminsearch(@FuncGrad_Aper_Synth_1D, c);
c = cnew;

%%%%%%%%%%%%%%%%
% FINAL RESULT %
%%%%%%%%%%%%%%%%
% --- Chebyshev
% f0 = InterpChebyshev(c);                                        % --- Coefficients of f0(xsi)
% f = expansiveFuncGen(f0);
% x = dmin / delta_xsi * polyval(f, xsi);
% x = x - min(x);

% --- Lagrange 
f0  = mapLagrangeInterpFunction(xsi_nodes, c);                    % --- Coefficients of f0(xsi)
f   = expansiveFuncGen(f0);
x   = polyval(f, xsi);

Spectrum_final = NFFT2_1D(Ea, x /(lambda / 2), Nu, 2, 6);
% Spectrum_final=Spectrum_final/max(abs(Spectrum_final));
   
handle = 0;

handle = handle + 1;
figure(handle)
plot(uu, 20 * log10(abs(Spectrum_final)   ./ max(abs(Spectrum_final))),   '-')
hold on
plot(uu, 20 * log10(abs(Spectrum_initial) ./ max(abs(Spectrum_initial))), 'm-.')
plot(uu, 20 * log10(abs(Spectrum)         ./ max(abs(Spectrum))),         'k-')
plot(uu, 20 * log10(abs(Upper_mask)),                                     'g')
plot(uu, 20 * log10(abs(Lower_mask)),                                     'r')
% axis([-pi pi -60 20]), xlabel('u'), ylabel('Spectrum_final (dB)')
hold off
drawnow

handle = handle + 1;
figure(handle)
plot(x / lambda, 1.2*ones(size(x)), 'x-', x_tapered_aperture / lambda, ones(size(x_tapered_aperture)), 'rx-', xstart / lambda, .8*ones(size(xstart)), 'k-x')
axis([min([min(x / lambda), min(xstart / lambda), min(x_tapered_aperture / lambda)]), max([max(x / lambda), max(xstart / lambda), max(x_tapered_aperture / lambda)]), .5, 1.5]);  
legend('final', 'tapered aperture', 'initial', 'Location', 'NorthWest')
xlabel('x / \lambda')

save Intermediate_result.mat c

