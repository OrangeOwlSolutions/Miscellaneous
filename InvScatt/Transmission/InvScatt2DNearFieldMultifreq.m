clear all
close all
clc
warning off

global leg leg_der N_sampling_points a NDF_x NDF_z d X_prime Z_prime prolate_spheroidals_x prolate_spheroidals_z beta numFrequencies transmission_reflection
 
c       = 3e8;                                                      % --- Light speed
   
lambda0 = 1;                                                        % --- Wavelength at the center-band

a_prime = 2 * lambda0;                                              % --- Half-dimension of the investigation domain along x
b_prime = 2 * lambda0;                                              % --- Half-dimension of the investigation domain along z

% numFrequencies = 15;                                                % --- Number of frequencies
numFrequencies = 15;                                                % --- Number of frequencies
Deltaf  = c / (4 * b_prime);                                        % --- Sampling step in frequency
freq    = (-(numFrequencies - 1) / 2 : (numFrequencies - 1) / 2) * Deltaf + c / lambda0;
                                                                    % --- Frequencies

lambda  = c ./ freq;                                                % --- Wavelengths
beta    = 2 * pi ./ lambda;                                         % --- Wavenumbers

d       = 4 * max(lambda);                                          % --- Distance between the investigation domain and the measurement domain
a       = 7 * max(lambda);                                          % --- Half-dimension of the measurement domain

SNR     = 30;                                                       % --- Signal to Noise Ratio                 

transmission_reflection = -1;                                        % --- -1 for transmission, 1 for reflection

handle = 0;                                                         % --- Pre-initial figure handle

d_d = gpuArray(d);
d_beta = gpuArray(beta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCRETIZATION OF THE INVESTIGATION DOMAIN %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Delta_x_prime       = min(lambda) / 6;                              % --- Discretization sampling step along x
Delta_z_prime       = min(lambda) / 6;                              % --- Discretization sampling step along z

Ppr_x               = ceil(a_prime / Delta_x_prime);                % --- 2 * Ppr_x + 1 discretization points along x
Ppr_z               = ceil(b_prime / Delta_z_prime);                % --- 2 * Ppr_z + 1 discretization points along z

x_prime             = a_prime * linspace(-1, 1, 2 * Ppr_x + 1);     % --- Sampling points along x
z_prime             = b_prime * linspace(-1, 1, 2 * Ppr_z + 1);     % --- Sampling points along z
[X_prime, Z_prime]  = meshgrid(x_prime, z_prime);                   % --- Sampling grid

d_X_prime = gpuArray(X_prime);
d_Z_prime = gpuArray(Z_prime);

%%%%%%%%%%%%%%%%%%%
% PSWF GENERATION %
%%%%%%%%%%%%%%%%%%%

NDF_x = ceil(1   * 4 * a_prime / min(lambda));                      % --- Number of prolate spheroidals along x
NDF_z = ceil(1   * 4 * b_prime / min(lambda));                      % --- Number of prolate spheroidals along z

prolate_spheroidals_x = Gen_prolate(a_prime, 2 * Ppr_x + 1);        % --- Generate prolate spheroidal wavefunctions along x at 2 * Ppr_x + 1 points
prolate_spheroidals_z = Gen_prolate(b_prime, 2 * Ppr_z + 1);        % --- Generate prolate spheroidal wavefunctions along z at 2 * Ppr_z + 1 points

prolate_spheroidals_x = prolate_spheroidals_x(1 : NDF_x, :);        % --- Take the first NDF_x PSWFs along x
prolate_spheroidals_z = prolate_spheroidals_z(1 : NDF_z, :);        % --- Take the first NDF_z PSWFs along z

%%%%%%%%%%%%%%%%%%%
% REFERENCE IMAGE %
%%%%%%%%%%%%%%%%%%%

Reference_image = zeros(size(X_prime));
Reference_image(find((abs(X_prime) < 0.15 * a_prime) & (abs(Z_prime) < 0.15 * b_prime))) = 1; 

handle = handle + 1;
figure(handle)
imagesc(z_prime, x_prime, abs(Reference_image.')), colorbar
xlabel('z [m]'), ylabel('x [m]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RADIATION OPERATOR: UNIFORM SAMPLING %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_uniform = -a : min(lambda) / 2: a;
N_sampling_points_uniform = length(x_uniform);

Radiation_Operator_Uniform = zeros(N_sampling_points_uniform * numFrequencies, NDF_x * NDF_z);

for s = 1 : numFrequencies,
    for k = 1 : N_sampling_points_uniform,
    
        R = sqrt((x_uniform(k) - X_prime).^2 + (d - Z_prime).^2);
   
        for p = 1 : NDF_x,
            for q = 1 : NDF_z,
                [PROLATE_X, PROLATE_Z] = meshgrid(prolate_spheroidals_x(p, :), prolate_spheroidals_z(q, :).');
                PROLATE_XZ = PROLATE_X .* PROLATE_Z;
                Radiation_Operator_Uniform((s - 1) * (N_sampling_points_uniform) + k, (p - 1) * NDF_z + q) = ...
                               sum(sum((besselh(0, 2, beta(s) * R) .* exp(transmission_reflection * 1i * beta(s) * Z_prime) .* PROLATE_XZ)));
            end
        end 
        100 * ((s - 1) * N_sampling_points_uniform + k) / (N_sampling_points_uniform * numFrequencies)

    end
end

S_uniform = svd(Radiation_Operator_Uniform);

handle = handle + 1;
plot(20 * log10(S_uniform / S_uniform(1)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECONSTRUCTION: UNIFORM SAMPLING %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data_uniform = zeros(size(Radiation_Operator_Uniform, 1), 1);
for s = 1 : numFrequencies,
    for k = 1 : N_sampling_points_uniform,
        R = sqrt((x_uniform(k) - X_prime).^2 + (d - Z_prime).^2);
        Data_uniform((s - 1) * (N_sampling_points_uniform) + k) = sum(sum((besselh(0, 2, beta(s) * R) .* exp(transmission_reflection * 1i * beta(s) * Z_prime) .* Reference_image)));
    end
end

Noise = randn(size(Data_uniform)) + 1i * randn(size(Data_uniform));
Noise = 10^(-SNR / 20) * sqrt(sum(sum(abs(Data_uniform).^2))) * Noise / sqrt(sum(sum(abs(Noise).^2)));
Data_uniform = Data_uniform + Noise;

Reconstruction_uniform = zeros(size(PROLATE_XZ));
Reconstruction_coefficients_uniform = pinv(Radiation_Operator_Uniform, 10^(-SNR / 20) * S_uniform(1)) * Data_uniform;
Reconstruction_coefficients_uniform = reshape(Reconstruction_coefficients_uniform, NDF_x, NDF_z);
for p = 1 : NDF_x,
    for q = 1 : NDF_z,
        [PROLATE_X, PROLATE_Z] = meshgrid(prolate_spheroidals_x(p, :), prolate_spheroidals_z(q, :).');
        PROLATE_XZ = PROLATE_X .* PROLATE_Z;
        Reconstruction_uniform = Reconstruction_uniform + Reconstruction_coefficients_uniform(q, p) * PROLATE_XZ;
    end
end

handle = handle + 1;
figure(handle)
imagesc(z_prime, x_prime, abs(Reconstruction_uniform.') / max(max(abs(Reconstruction_uniform)))), colorbar
xlabel('z [m]'), ylabel('x [m]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STARTING SAMPLING POINTS FOR THE OPTIMIZATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_unknowns = 3;                                                             % --- Number of Legendre polynomials

N_min_sampling_points = 50;
N_max_sampling_points = 50;

Coeff_legendre_opt      = zeros(1, N_unknowns);
% Coeff_legendre_opt(1)   = (min(lambda) / 2) * ((N_min_sampling_points - 1) / 2) / a;    % --- Initial uniform spacing equal to lambda_min / 2
Coeff_legendre_opt(1)   = (lambda0 / 2) * ((N_min_sampling_points - 1) / 2) / a;    % --- Initial uniform spacing equal to lambda_min / 2

% load Coeff_temp_opt.mat
load Coeff_opt_final_50_points.mat

Func_values = zeros(1, N_max_sampling_points);

for N_sampling_points = N_min_sampling_points : 2: N_max_sampling_points,

    % --- Forcing the number of sampling points to be even
    if mod(N_sampling_points, 2) ~= 0
        N_sampling_points = N_sampling_points + 1;                              
    end

    x_normalized = linspace(0, 1, N_sampling_points / 2);
    % --- Constructing the basis functions for the point represenation
    leg     = zeros(N_unknowns, length(x_normalized));
    leg_der = zeros(N_unknowns, length(x_normalized));
    for k = 1 : N_unknowns,
        % --- Legendre polynomials
        leg(k, :)       = polyval(LegendrePoly(k),          x_normalized);
    end
    for k = 1 : N_unknowns,
        % --- Legendre polynomial derivatives
        leg_der(k, :)   = polyval(polyder(LegendrePoly(k)),          x_normalized);
    end

    x_start_optimization        = zeros(1, N_sampling_points / 2);
    x_start_optimization_der    = zeros(1, N_sampling_points / 2);
    for k = 2 : N_unknowns,
    end
    for k = 1 : N_unknowns,
        x_start_optimization        = x_start_optimization      + Coeff_legendre_opt(k) * leg(k, :);
        x_start_optimization_der    = x_start_optimization_der  + Coeff_legendre_opt(k) * leg_der(k, :);
    end
    % --- Enforcing symmetry around zero
    if (min(x_start_optimization_der) < 0)
        x_start_optimization = x_start_optimization - min(x_start_optimization_der) * leg(1, :);
    end
    temp = x_start_optimization;
    temp = temp - min(temp);
    temp = temp(2 : length(temp));
    x_start_optimization = [-fliplr(temp) 0 temp];
    x_start_optimization = a * x_start_optimization / max(x_start_optimization);

    handle = handle + 1;
    figure(handle)
    plot(x_start_optimization, zeros(size(x_start_optimization)), 'o')
    hold on

    %%%%%%%%%%%%%%%%%%%%%%
    % INITIAL FUNCTIONAL %
    %%%%%%%%%%%%%%%%%%%%%%

    Radiation_Operator_Init = zeros(N_sampling_points * numFrequencies, NDF_x * NDF_z);

    for s = 1 : numFrequencies,
        for k = 1 : N_sampling_points - 1,
    
            R = sqrt((x_start_optimization(k) - X_prime).^2 + (d - Z_prime).^2);
   
            for p = 1 : NDF_x,
                for q = 1 : NDF_z,
                    [PROLATE_X, PROLATE_Z] = meshgrid(prolate_spheroidals_x(p, :), prolate_spheroidals_z(q, :).');
                    PROLATE_XZ = PROLATE_X .* PROLATE_Z;
                    Radiation_Operator_Init((s - 1) * (N_sampling_points - 1) + k, (p - 1) * NDF_z + q) = sum(sum((besselh(0, 2, beta(s) * R) .* exp(transmission_reflection * 1i * beta(s) * Z_prime) .* PROLATE_XZ)));
                end
            end 
            100 * ((s - 1) * N_sampling_points + k) / (N_sampling_points * numFrequencies)
        end
    end

    S_init = svd(Radiation_Operator_Init);

    handle = handle + 1;
    handle_singular_values = handle;
    figure(handle)
    plot(20 * log10(S_init / S_init(1)))

    %%%%%%%%%%%%%%%%
    % OPTIMIZATION %
    %%%%%%%%%%%%%%%%

%     [Coeff_legendre_opt, Func_values(N_sampling_points)] = fminunc(@Inverse_Scattering_2D_Near_Field_Multifrequency_main_Functional, Coeff_legendre_opt);

end

%%%%%%%%%%%%%%%%%
% FINAL RESULTS %
%%%%%%%%%%%%%%%%%

x_opt       = zeros(1, N_sampling_points / 2);
x_opt_der   = zeros(1, N_sampling_points / 2);
for k = 1 : N_unknowns,
    x_opt       = x_opt     + Coeff_legendre_opt(k) * leg(k, :);
    x_opt_der   = x_opt_der + Coeff_legendre_opt(k) * leg_der(k, :);
end
if (min(x_opt_der) < 0)
    x_opt = x_opt - min(x_opt_der) * leg(1, :);
end
temp = x_opt;
temp = temp - min(temp);
temp = temp(2 : length(temp));
x_opt = [-fliplr(temp) 0 temp];
x_opt = a * x_opt / max(abs(x_opt));

N_sampling_points = N_sampling_points - 1;
Radiation_Operator_Opt = zeros(N_sampling_points * numFrequencies, NDF_x * NDF_z);

for s = 1 : numFrequencies,
    for k = 1 : N_sampling_points,
    
        R = sqrt((x_opt(k) - X_prime).^2 + (d - Z_prime).^2);
   
        for p = 1 : NDF_x,
            for q = 1 : NDF_z,
                [PROLATE_X, PROLATE_Z] = meshgrid(prolate_spheroidals_x(p, :), prolate_spheroidals_z(q, :).');
                PROLATE_XZ = PROLATE_X .* PROLATE_Z;
                Radiation_Operator_Opt((s - 1) * (N_sampling_points) + k, (p - 1) * NDF_z + q)=sum(sum((besselh(0, 2, beta(s) * R) .* exp(j * beta(s) * Z_prime) .* PROLATE_XZ)));
            end
        end 
        100 * ((s - 1) * N_sampling_points + k) / (N_sampling_points * numFrequencies)
    end
end

S_opt = svd(Radiation_Operator_Opt);

handle = handle + 1;
figure(handle)
plot(x_opt, zeros(size(x_opt)), 'o')
hold on

handle = handle + 1;
handle_singular_values = handle;
figure(handle)
plot(20 * log10(S_opt / S_opt(1)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECONSTRUCTION: NON-UNIFORM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data_non_uniform = zeros(size(Radiation_Operator_Opt, 1), 1);
for s = 1 : numFrequencies,
    for k = 1 : N_sampling_points,
        R = sqrt((x_opt(k) - X_prime).^2 + (d - Z_prime).^2);
        Data_non_uniform((s - 1) * (N_sampling_points) + k) = sum(sum((besselh(0, 2, beta(s) * R) .* exp(transmission_reflection * 1i * beta(s) * Z_prime) .* Reference_image)));
    end
end

Noise = randn(size(Data_non_uniform)) + 1i * randn(size(Data_non_uniform));
Noise = 10^(-SNR / 20) * sqrt(sum(sum(abs(Data_non_uniform).^2))) * Noise / sqrt(sum(sum(abs(Noise).^2)));
Data_non_uniform = Data_non_uniform + Noise;

Reconstruction_non_uniform = zeros(size(PROLATE_XZ));
Reconstruction_coefficients_non_uniform = pinv(Radiation_Operator_Opt, 10^(-SNR / 20) * S_opt(1)) * Data_non_uniform;
Reconstruction_coefficients_non_uniform = reshape(Reconstruction_coefficients_non_uniform, NDF_x, NDF_z);
for p = 1 : NDF_x,
    for q = 1 : NDF_z,
        [PROLATE_X, PROLATE_Z] = meshgrid(prolate_spheroidals_x(p, :), prolate_spheroidals_z(q, :).');
        PROLATE_XZ = PROLATE_X .* PROLATE_Z;
        Reconstruction_non_uniform = Reconstruction_non_uniform + Reconstruction_coefficients_non_uniform(q, p) * PROLATE_XZ;
    end
end

handle = handle + 1;
figure(handle)
imagesc(z_prime, x_prime, abs(Reconstruction_non_uniform.') / max(max(abs(Reconstruction_non_uniform)))), colorbar
xlabel('z [m]'), ylabel('x [m]')

save Coeff_temp_opt.mat Coeff_legendre_opt
