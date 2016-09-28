close all
clear all
clc
warning off

lambda  = 1;                                                        % --- Wavelength
beta    = 2 * pi / lambda;                                          % --- Wavenumber

d = 7 * lambda;                                                     % --- Distance of the measurement plane

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AP. PLANE DISCRETIZATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aap = 7 * lambda;                                                   % --- Half-size of the ap plane along x
bap = 7 * lambda;                                                   % --- Half-size of the ap plane along y

x_interp = 0 : lambda / 5 : aap;                                    % --- Discretization points of the ap plane along x
x_interp = [-fliplr(x_interp(2 : length(x_interp))) x_interp];
y_interp = 0 : lambda / 5 : bap;                                    % --- Discretization points of the ap plane along y
y_interp = [-fliplr(y_interp(2 : length(y_interp))) y_interp];
[X_interp, Y_interp] = meshgrid(x_interp, y_interp);

Nx_interp = length(x_interp);                                       % --- Number of discretization points of the ap plane along x
Ny_interp = length(y_interp);                                       % --- Number of discretization points of the ap plane along x
deltaX_interp = 2 * aap / Nx_interp;                                % --- Sampling step of the ap plane discretization along x
deltaY_interp = 2 * bap / Ny_interp;                                % --- Sampling step of the ap plane discretization along y
delta_x_interp = deltaX_interp * ones(1, Nx_interp);
delta_y_interp = deltaY_interp * ones(1, Ny_interp);
 
apSupport = ones(size(X_interp));                                   % --- Ap support
indices1 = find(sqrt(X_interp.^2 + Y_interp.^2) > aap);
apSupport(indices1) = zeros(size(indices1));     

c = beta * aap;                                                     % --- SBP

Nmax = 2;
nmax = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEASUR. PLANE DISCRETIZATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 10 * lambda;                                                    % --- Half-size of the measur plane along x
b = 10 * lambda;                                                    % --- Half-size of the measur plane along y

x = 0 : lambda / 2 : a;                                             % --- Discretization points of the measur plane along x
x = [-fliplr(x(2 : length(x))) x];
y = 0 : lambda / 2 : b;                                             % --- Discretization points of the measur plane along x
y = [-fliplr(y(2 : length(y))) y];
[X, Y] = meshgrid(x, y);

Nx = length(x);                                                     % --- Number of discretization points of the measur plane along x
Ny = length(y);                                                     % --- Number of discretization points of the measur plane along y
deltax = 2 * a / Nx;                                                % --- Sampling step of the measur plane discretization along x
deltay = 2 * b / Ny;                                                % --- Sampling step of the measur plane discretization along y
delta_x_1 = deltax * ones(1, Nx);
delta_y_1 = deltay * ones(1, Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRUM DISCRETIZATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

chi_uu = 6;                                                         % --- Oversampling factor along u
chi_vv = 6;                                                         % --- Oversampling factor along v

amax = max([aap a]);
bmax = max([bap b]);

uu = 0 : pi / (chi_uu * amax) : beta;                               % --- Spectrum discretization points along u
uu = [-fliplr(uu(2 : length(uu))) uu];
vv = 0 : pi / (chi_vv * bmax) : beta;                               % --- Spectrum discretization points along v
vv = [-fliplr(vv(2 : length(vv))) vv];
[UU, VV] = meshgrid(uu, vv);

Nu = length(uu);                                                    % --- Number of spectrum discretization points along u
Nv = length(vv);                                                    % --- Number of spectrum discretization points along v
delta_uu = diff(uu);                                                % --- Spectrum sampling steps along u
delta_uu = [delta_uu, delta_uu(length(delta_uu))];                 
delta_vv = diff(vv);                                                % --- Spectrum sampling steps along v
delta_vv = [delta_vv, delta_vv(length(delta_vv))];

%%%%%%%%%%%%%%%%
% PROP. FACTOR %
%%%%%%%%%%%%%%%%
ArgumentSquareRoot = beta^2 - UU.^2 - VV.^2;
WW = sqrt(beta^2 - UU.^2 - VV.^2);
WW(ArgumentSquareRoot < 0) = -1i * sqrt(UU(ArgumentSquareRoot < 0).^2 + VV(ArgumentSquareRoot < 0).^2 - beta^2);
Filter = zeros(size(WW));
indices = find(sqrt(UU.^2 + VV.^2) <= 0.99 * beta);
Filter(indices) = ones(size(indices)); 
PropFactor = Filter .* exp(-1i * d * WW);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINIZIONE CAMPO DI APERTURA, SPETTRO E CAMPO SUL PIANO DI MISURA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Reference ap field
EaRef = apSupport .* exp(-(20 / (20 * aap^2 * log10(exp(1)))) * (X_interp.^2 + Y_interp.^2));

% --- Reference spectrum
SpectrumRef = Filter .* Ea2S(uu, vv, EaRef, x_interp, y_interp, delta_x_interp, delta_y_interp);

% --- Measur plane field
V       = S2Ea(uu, vv, SpectrumRef .* PropFactor, x, y, delta_uu, delta_vv);                   
SNR     = 20;
Noise   = randn(size(V)) + 1i * randn(size(V));
Noise   = 10^(-SNR / 20) * sqrt(sum(sum(abs(V).^2))) * Noise / sqrt(sum(sum(abs(Noise).^2)));
V       = V + Noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIAL FUNCTIONS SECTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NDF_x   = ceil(4 * aap / lambda);
NDF_y   = ceil(4 * bap / lambda);

x_ap    = linspace(-aap, aap, 2 * NDF_x + 1);                       % --- Discretization points of the ap plane along x (for inversion)
y_ap    = linspace(-bap, bap, 2 * NDF_y + 1);                       % --- Discretization points of the ap plane along y (for inversion)
[X_ap, Y_ap] = meshgrid(x_ap, y_ap);

apSupport_ap = ones(size(X_ap));                                    % --- Ap support (for inversion)
indices2     = find(sqrt(X_ap.^2 + Y_ap.^2) > aap);
apSupport_ap(indices2) = zeros(size(indices2));     

PSI = computePSINn(c, Nmax, nmax, x_ap / aap, y_ap / bap);
PSI_interp = computePSINn(c, Nmax, nmax, x_interp / aap, y_interp / bap);

specialFunctions = zeros((2 * NDF_x + 1) * (2 * NDF_y + 1), Nmax * nmax);
specialFunctions_interp=zeros(length(x_interp) * length(y_interp), Nmax * nmax);

for p = 1 : Nmax,
    for q = 1 : nmax,
        specialFunctions(:, (p - 1) * nmax + q) = reshape(squeeze(PSI(p, q, :, :)),(2 * NDF_x + 1) * (2 * NDF_y + 1), 1);
        specialFunctions_interp(:, (p - 1) * nmax + q) = reshape(squeeze(PSI_interp(p, q, :, :)), length(x_interp) * length(y_interp), 1);
    end
end
        
SNR_SVD                 = 50;
S                       = svd(specialFunctions);
SVD_Threshold           = S(1) * 10^(-SNR_SVD / 20);
InvSpecialFunctions     = pinv(specialFunctions, SVD_Threshold);

clear specialFunctions

%%%%%%%%%%%%%%%%%%
% TRANSFORMATION %
%%%%%%%%%%%%%%%%%%

SpectrumNoFiltering         = Ea2S(uu, vv, V, x, y, delta_x_1, delta_y_1) .* conj(PropFactor);

EaFiltering                 = apSupport_ap .* S2Ea(uu, vv, SpectrumNoFiltering, x_ap, y_ap, delta_uu, delta_vv);        
EaFiltering_interp          = apSupport .* S2Ea(uu, vv, SpectrumNoFiltering, x_interp, y_interp, delta_uu, delta_vv);   

SpectrumFiltering           = Filter .* Ea2S(uu, vv, EaFiltering_interp, x_interp, y_interp, delta_x_interp, delta_y_interp);

C                           = InvSpecialFunctions * reshape(EaFiltering, numel(EaFiltering), 1);

EaSpecialFunctions          = reshape(specialFunctions_interp * C, length(x_interp), length(y_interp));
SpectrumSpecialFunctions    = Filter .* Ea2S(uu, vv, EaSpecialFunctions, x_interp, y_interp, delta_x_interp, delta_y_interp);

SpectrumRef                 = WW .* SpectrumRef;
SpectrumSpecialFunctions    = WW .* SpectrumSpecialFunctions;
SpectrumFiltering           = WW .* SpectrumFiltering;
SpectrumNoFiltering         = WW .* SpectrumNoFiltering;

%%%%%%%%%%
% GRAPHS %
%%%%%%%%%%

handle = 0;

% --- Spectrum
handle = handle + 1;
figure(handle)
plot(uu / beta, 20 * log10(abs(SpectrumRef((Nu - 1) / 2 + 1, :)) / max(abs(SpectrumRef((Nu - 1) / 2 + 1, :)))), '-.', ...
     uu / beta, 20 * log10(abs(SpectrumSpecialFunctions((Nu - 1) / 2 + 1,:)) / max(abs(SpectrumSpecialFunctions((Nu - 1) / 2 + 1, :)))), 'r', ...
     uu / beta, 20 * log10(abs(SpectrumFiltering((Nu - 1) / 2 + 1,:)) / max(abs(SpectrumFiltering((Nu - 1) / 2 + 1, :)))), 'k', 'LineWidth', 2)
xlabel('k_x/\beta'), ylabel('Far-Field amplitude [dB]')
axis([-1 1 -80 0])

handle = handle + 1;
figure(handle)
plot(vv / beta, 20 * log10(abs(SpectrumRef(:, (Nv - 1) / 2 + 1)) / max(abs(SpectrumRef(:, (Nv - 1) / 2 + 1)))), '-.', ...
     vv / beta, 20 * log10(abs(SpectrumSpecialFunctions(:, (Nv - 1) / 2 + 1)) / max(abs(SpectrumSpecialFunctions(:, (Nv - 1) / 2 + 1)))), 'r', ...
     vv / beta, 20 * log10(abs(SpectrumFiltering(:, (Nv - 1) / 2 + 1)) / max(abs(SpectrumFiltering(:, (Nv - 1) / 2 + 1)))), 'k', 'LineWidth', 2)
xlabel('k_y/\beta'), ylabel('Far-Field amplitude [dB]')
axis([-1 1 -80 0])

handle = handle + 1;
figure(handle)
contourf(uu / beta, uu / beta, 20 * log10(abs(SpectrumRef) / max(max(abs(SpectrumRef)))),[-70 -60 -50 -40 -30 -20 -15 -10 -5 -7 -3 -1 0])
xlabel('k_x/\lambda')
ylabel('k_y/\lambda')
colorbar

handle = handle + 1;
figure(handle)
contourf(uu / beta, uu / beta, 20 * log10(abs(SpectrumSpecialFunctions) / max(max(abs(SpectrumSpecialFunctions)))),[-70 -60 -50 -40 -30 -20 -15 -10 -5 -7 -3 -1 0])
xlabel('k_x/\lambda')
ylabel('k_y/\lambda')
colorbar

handle = handle + 1;
figure(handle)
contourf(uu / beta, uu / beta, 20 * log10(abs(SpectrumFiltering) / max(max(abs(SpectrumFiltering)))),[-70 -60 -50 -40 -30 -20 -15 -10 -5 -7 -3 -1 0])
xlabel('k_x/\lambda')
ylabel('k_y/\lambda')
colorbar

handle = handle + 1;
figure(handle)
contourf(uu / beta, uu / beta, 20 * log10(abs(SpectrumNoFiltering) / max(max(abs(SpectrumNoFiltering)))),[-70 -60 -50 -40 -30 -20 -15 -10 -5 -7 -3 -1 0])
xlabel('k_x/\lambda')
ylabel('k_y/\lambda')
colorbar

% --- Ap field
handle = handle + 1;
figure(handle)
contourf(x_interp / lambda, y_interp / lambda, 20 * log10(abs(EaRef) / max(max(abs(EaRef)))),[-20 -15 -10 -5 -7 -3 -1 0])
xlabel('x/\lambda')
ylabel('y/\lambda')
colorbar

handle = handle + 1;
figure(handle)
contourf(x_interp / lambda, y_interp / lambda, 20 * log10(abs(EaSpecialFunctions) / max(max(abs(EaSpecialFunctions)))),[-20 -15 -10 -5 -7 -3 -1 0])
xlabel('x/\lambda')
ylabel('y/\lambda')
colorbar

handle = handle + 1;
figure(handle)
contourf(x_interp / lambda, y_interp / lambda, 20 * log10(abs(EaFiltering_interp) / max(max(abs(EaFiltering_interp)))),[-20 -15 -10 -5 -7 -3 -1 0])
xlabel('x/\lambda')
ylabel('y/\lambda')
colorbar

