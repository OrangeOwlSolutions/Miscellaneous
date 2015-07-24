function data = FBP(data, a, b)

global c                                                        % --- Speed of light [m/s]
global fc                                                       % --- Center frequency [Hz]
global BW                                                       % --- Bandwidth [Hz]
global f1                                                       % --- Minimum frequency [Hz]
global f2                                                       % --- Maximum frequency [Hz]

OverSamplingFactorRho       = 2;                                % --- Oversampling factor in the rho variable
OverSamplingFactorCosAlpha  = 4;                                % --- Oversampling factor in the cosalpha variable

Nextra = 32;                                                    % --- Number of extra samples to prevent errors at the borders

ExtrapVal = nan;

% ??? Cambiare secondo Yegulalp
N_PulsesPerSubAperture = 2;
N_SubApertures  = ((b - a + 1) / N_PulsesPerSubAperture);

if(mod((b - a + 1), N_PulsesPerSubAperture) ~= 0)
    error('The aperture length must be an integer number of Sub-Aperture length')
end

% --- Partial image initialization with all zeros
data.Nearest    = zeros(size(data.x_mat));
data.Linear     = zeros(size(data.x_mat));
data.Cubic      = zeros(size(data.x_mat));
data.Spline     = zeros(size(data.x_mat));
data.knabNFFT   = zeros(size(data.x_mat));
data.rectNFFT   = zeros(size(data.x_mat));
data.knab       = zeros(size(data.x_mat));
    
% --- Loop through every subaperture
for s = 0 : N_SubApertures - 1
    
    fprintf ('Sub-aperture %d of %d\n', s+1, N_SubApertures);
    
    % --- Subaperture center
    center = s * N_PulsesPerSubAperture + a;
    
    % --- Creates the polar grid
    % --- \rho_m  = |\underline{r}  - \underline{r}_a(t_m)|
    RHO_OUT      = sqrt((data.x_mat - data.AntX(center)).^2 + (data.y_mat - data.AntY(center)).^2 + (data.z_mat - data.AntZ(center)).^2);
    % --- \cos \alpha_m = \frac{\left(\underline{r} - \underline{r}_a(t_m)\right)\cdot \left(\underline{r}_a(\xi_n+t_m)  - \underline{r}_a(t_m)\right)}{\rho_m\xi_n}
    Dxs = (data.AntX(center + 1) - data.AntX(center));
    Dys = (data.AntY(center + 1) - data.AntY(center));
    Dzs = (data.AntZ(center + 1) - data.AntZ(center));
    dxs = Dxs / sqrt(Dxs^2 + Dys^2 + Dzs^2);
    dys = Dys / sqrt(Dxs^2 + Dys^2 + Dzs^2);
    dzs = Dzs / sqrt(Dxs^2 + Dys^2 + Dzs^2);
    COSALPHA_OUT = ((data.x_mat - data.AntX(center)) ./ RHO_OUT * dxs + ...
                    (data.y_mat - data.AntY(center)) ./ RHO_OUT * dys + ...
                    (data.z_mat - data.AntZ(center)) ./ RHO_OUT * dzs);
    
    % --- Maximum and minimum values of rho_out and cos(alpha)_out
    maxrho      = max(max(RHO_OUT));
    minrho      = min(min(RHO_OUT));
    maxcosalpha = max(max(COSALPHA_OUT));
    mincosalpha = min(min(COSALPHA_OUT));

    SubApertureLength = sqrt((data.AntX(center + N_PulsesPerSubAperture - 1) - data.AntX(center)).^2 + ...
                             (data.AntY(center + N_PulsesPerSubAperture - 1) - data.AntY(center)).^2 + ...
                             (data.AntZ(center + N_PulsesPerSubAperture - 1) - data.AntZ(center)).^2);
    
    % --- Number of samples in the output rho_out and cos(alpha)_out variables 
    Nrho        = (floor((maxrho - minrho) / (c / (2 * BW))) + 1) * 2 * OverSamplingFactorRho;
    Ncosalpha   = (floor((maxcosalpha - mincosalpha) / (c / (2 * f2 * SubApertureLength))) + 1) * 2 * OverSamplingFactorCosAlpha;
    
    % --- Sampling steps in the output rho_in and alpha variables 
    deltarho      = (maxrho-minrho)/(Nrho-1);
    deltacosalpha = (maxcosalpha-mincosalpha)/(Ncosalpha-1);
    
    % --- Output rho_in and alpha variables
    rho_in       = (-Nextra: Nrho   + Nextra - 1) * deltarho      + minrho;
    cosalpha_in  = (-Nextra: Ncosalpha + Nextra - 1) * deltacosalpha + mincosalpha;
    
    [RHO_IN, COSALPHA_IN] = meshgrid(rho_in, cosalpha_in);
    
    % --- BackProjection on the output local polar grid
    data_temp               = data;
    data_temp.pho_mat       = RHO_IN;
    data_temp.COSALPHA_OUT  = COSALPHA_IN;
    data_temp = BP_POLAR_NFFT(data_temp, center, center + N_PulsesPerSubAperture - 1);
    data_temp.im_final = data_temp.im_final .* exp(-1j * 2 * pi * (2 * (RHO_IN - data.R0(a)) / c) * fc);

    % --- Interpolation of the image from the input local polar grid to the final global polar grid
    data.Nearest  = data.Nearest  + interp2(RHO_IN, COSALPHA_IN, data_temp.im_final, RHO_OUT, COSALPHA_OUT, 'nearest', ExtrapVal) .* ...
                                    exp(1j * 2 * pi * 2 * (RHO_OUT - data.R0(a)) / c * fc);
    data.Linear   = data.Linear   + interp2(RHO_IN, COSALPHA_IN, data_temp.im_final, RHO_OUT, COSALPHA_OUT, 'linear', ExtrapVal) .* ...
                                    exp(1j * 2 * pi * 2 * (RHO_OUT - data.R0(a)) / c * fc);
    data.Cubic    = data.Cubic    + interp2(RHO_IN, COSALPHA_IN, data_temp.im_final, RHO_OUT, COSALPHA_OUT, 'cubic', ExtrapVal) .* ...
                                    exp(1j * 2 * pi * 2 * (RHO_OUT - data.R0(a)) / c * fc);
    data.Spline   = data.Spline   + interp2(RHO_IN, COSALPHA_IN, data_temp.im_final, RHO_OUT, COSALPHA_OUT, 'spline', ExtrapVal) .* ...
                                    exp(1j * 2 * pi * 2 * (RHO_OUT - data.R0(a)) / c * fc);
    data.knabNFFT = data.knabNFFT + knabNFFT2D(ExtrapVal, RHO_IN, COSALPHA_IN, data_temp.im_final, RHO_OUT, COSALPHA_OUT, (1 / deltarho) * 0.5, (1 / deltacosalpha) * 0.5) .* ...
                                    exp(1j * 2 * pi * 2 * (RHO_OUT - data.R0(a)) / c * fc);
    data.rectNFFT = data.rectNFFT + rectNFFT2D(ExtrapVal, RHO_IN, COSALPHA_IN, data_temp.im_final, RHO_OUT, COSALPHA_OUT, (1 / deltarho) * 0.5, (1 / deltacosalpha) * 0.5) .* ...
                                    exp(1j * 2 * pi * 2 * (RHO_OUT - data.R0(a)) / c * fc);
    data.knab     = data.knab     + knab2D(ExtrapVal, RHO_IN, COSALPHA_IN, data_temp.im_final, RHO_OUT, COSALPHA_OUT, (1 / deltarho) * 0.5, (1 / deltacosalpha) * 0.5) .* ...
                                    exp(1j * 2 * pi * 2 * (RHO_OUT - data.R0(a)) / c * fc);
    
end



end

