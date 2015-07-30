function data = FFBP_POLAR_NFFT(data, a, b, level)

global c                                                        % --- Speed of light [m/s]

global fc                                                       % --- Center frequency [Hz]
global BW                                                       % --- Bandwidth [Hz]
global f1                                                       % --- Minimum frequency [Hz]
global f2                                                       % --- Maximum frequency [Hz]

ExtrapVal = nan;

fprintf('__________________\nlevel = %d\n____________________\n',level);

% --- Partial image initialization with all zeros
data.Nearest    = zeros(size(data.RHO_OUT));
data.Linear     = zeros(size(data.RHO_OUT));
data.Cubic      = zeros(size(data.RHO_OUT));
data.Spline     = zeros(size(data.RHO_OUT));
data.knabNFFT   = zeros(size(data.RHO_OUT));
data.rectNFFT   = zeros(size(data.RHO_OUT));
data.knab       = zeros(size(data.RHO_OUT));

% --- Create two sub-apertures
center = floor((a+b)/2);

% --- Short aperture: stop the recursion
if(b-a <= 2)
    data.pho_mat = data.RHO_OUT;
    data            = BP_POLAR_NFFT(data,a,b);
    data.Nearest    = data.im_final;
    data.Linear     = data.im_final;
    data.Cubic      = data.im_final;
    data.Spline     = data.im_final;
    data.knabNFFT   = data.im_final;
    data.rectNFFT   = data.im_final;
    data.knab       = data.im_final;
else
    
    SINALPHA_OUT = sqrt(1-data.COSALPHA_OUT.^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIRST APERTURE [a, center] %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dist1 =    sqrt((data.AntX(a) - data.AntX(a))^2 + (data.AntY(a) - data.AntY(a))^2 + (data.AntZ(a) - data.AntZ(a))^2);
    RHO_OUT1 = sqrt((data.RHO_OUT .* data.COSALPHA_OUT - dist1).^2 + (data.RHO_OUT .* SINALPHA_OUT).^2);
    COSALPHA_OUT1 = cos(atan2(data.RHO_OUT .* SINALPHA_OUT, data.RHO_OUT .* data.COSALPHA_OUT - dist1));

    % --- Maximum and minimum values of rho_out and cos(alpha)_out
    maxrho1 = max(max(RHO_OUT1));
    minrho1 = min(min(RHO_OUT1));
    maxcosalpha1 = max(max(COSALPHA_OUT1));
    mincosalpha1 = min(min(COSALPHA_OUT1));

    SubApertureLength1 =  sqrt((data.AntX(center) - data.AntX(a)).^2 + (data.AntY(center) - data.AntY(a)).^2 + (data.AntZ(center) - data.AntZ(a)).^2);

    % --- Number of samples in the output rho_out and cos(alpha)_out variables 
    % ???
    Nrho1      = floor((maxrho1 - minrho1) / (c / (2 * BW))) * 5;
    Ncosalpha1 = floor((maxcosalpha1 - mincosalpha1) / (c / (2 * f2 * SubApertureLength1))) * 5 + 3;

    % --- Create the new angular sparse image grid for the first sub aperture
    % ???
    rho_in1        = linspace(minrho1 - (30 * (maxrho1 - minrho1) / (Nrho1 - 1)), maxrho1 + (30 * (maxrho1 - minrho1) / (Nrho1 - 1)), Nrho1 +  60);
    cosalpha_in1   = linspace(mincosalpha1 - (30 * (maxcosalpha1 - mincosalpha1) / (Ncosalpha1 - 1)), maxcosalpha1 + (30 * (maxcosalpha1 - mincosalpha1) / (Ncosalpha1 - 1)), Ncosalpha1 + 60);

    [RHO_IN1, COSALPHA_IN1] = meshgrid(rho_in1,cosalpha_in1);
    data_temp1 = data;
    data_temp1.RHO_OUT = RHO_IN1;
    data_temp1.COSALPHA_OUT = COSALPHA_IN1;

    data_temp1 = FFBP_POLAR_NFFT(data_temp1, a, center, level+1);

    data_temp1.Nearest  = data_temp1.Nearest    .* exp(-1j * (4 * pi / c * RHO_IN1) * fc);
    data_temp1.Linear   = data_temp1.Linear     .* exp(-1j * (4 * pi / c * RHO_IN1) * fc);
    data_temp1.Cubic    = data_temp1.Cubic      .* exp(-1j * (4 * pi / c * RHO_IN1) * fc);
    data_temp1.Spline   = data_temp1.Spline     .* exp(-1j * (4 * pi / c * RHO_IN1) * fc);
    data_temp1.knabNFFT = data_temp1.knabNFFT   .* exp(-1j * (4 * pi / c * RHO_IN1) * fc);
    data_temp1.rectNFFT = data_temp1.rectNFFT   .* exp(-1j * (4 * pi / c * RHO_IN1) * fc);
    data_temp1.knab     = data_temp1.knab       .* exp(-1j * (4 * pi / c * RHO_IN1) * fc);

    deltarho1       = mean(diff(RHO_IN1(1, :)));
    deltacosalpha1  = mean(diff(COSALPHA_IN1(:, 1)));

    % --- Interpolation of the image from the input local polar grid to the final global polar grid
    Reconstruction_Nearest1     = interp2(              RHO_IN1, COSALPHA_IN1, data_temp1.Nearest,  RHO_OUT1, COSALPHA_OUT1, 'nearest', ExtrapVal) .* exp(1j * (4 * pi / c * RHO_OUT1) * fc);
    Reconstruction_Linear1      = interp2(              RHO_IN1, COSALPHA_IN1, data_temp1.Linear,   RHO_OUT1, COSALPHA_OUT1, 'linear',  ExtrapVal) .* exp(1j * (4 * pi / c * RHO_OUT1) * fc);
    Reconstruction_Cubic1       = interp2(              RHO_IN1, COSALPHA_IN1, data_temp1.Cubic,    RHO_OUT1, COSALPHA_OUT1, 'cubic',   ExtrapVal) .* exp(1j * (4 * pi / c * RHO_OUT1) * fc);
    Reconstruction_Spline1      = interp2(              RHO_IN1, COSALPHA_IN1, data_temp1.Spline,   RHO_OUT1, COSALPHA_OUT1, 'spline',  ExtrapVal) .* exp(1j * (4 * pi / c * RHO_OUT1) * fc);
    Reconstruction_KnabNFFT1    = knabNFFT2D(ExtrapVal, RHO_IN1, COSALPHA_IN1, data_temp1.knabNFFT, RHO_OUT1, COSALPHA_OUT1, (1 / deltarho1) * 0.5, (1 / deltacosalpha1) * 0.5) .* exp(1j * (4 * pi / c * RHO_OUT1) * fc);
    Reconstruction_RectNFFT1    = rectNFFT2D(ExtrapVal, RHO_IN1, COSALPHA_IN1, data_temp1.rectNFFT, RHO_OUT1, COSALPHA_OUT1, (1 / deltarho1) * 0.5, (1 / deltacosalpha1) * 0.5) .* exp(1j * (4 * pi / c * RHO_OUT1) * fc);
    Reconstruction_Knab1        = knab2D(    ExtrapVal, RHO_IN1, COSALPHA_IN1, data_temp1.knab,     RHO_OUT1, COSALPHA_OUT1, (1 / deltarho1) * 0.5, (1 / deltacosalpha1) * 0.5) .* exp(1j * (4 * pi / c * RHO_OUT1) * fc);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SECOND APERTURE [center + 1, b] %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dist2 =    sqrt((data.AntX(center + 1) - data.AntX(a))^2 + (data.AntY(center + 1) - data.AntY(a))^2 + (data.AntZ(center + 1) - data.AntZ(a))^2);
    RHO_OUT2 = sqrt((data.RHO_OUT .* data.COSALPHA_OUT - dist2).^2 + (data.RHO_OUT .* SINALPHA_OUT).^2);
    COSALPHA_OUT2 = cos(atan2(data.RHO_OUT .* SINALPHA_OUT, data.RHO_OUT .* data.COSALPHA_OUT - dist2));

    % --- Maximum and minimum values of rho_out and cos(alpha)_out
    maxrho2 = max(max(RHO_OUT2));
    minrho2 = min(min(RHO_OUT2));
    maxcosalpha2 = max(max(COSALPHA_OUT2));
    mincosalpha2 = min(min(COSALPHA_OUT2));

    SubApertureLength2 =  sqrt((data.AntX(b) - data.AntX(center)).^2 + (data.AntY(b) - data.AntY(center)).^2 + (data.AntZ(b) - data.AntZ(center)).^2);

    % --- Number of samples in the output rho_out and cos(alpha)_out variables 
    % ???
    Nrho2       = floor((maxrho2 - minrho2) / (c / (2 * BW))) * 5;
    Ncosalpha2  = floor((maxcosalpha2 - mincosalpha2) / (c / (2 * f2 * SubApertureLength2))) * 5 + 3;

    % --- Create the new angular sparse image grid for the first sub aperture
    % ???
    rho_in2         = linspace(minrho2 - (30 * (maxrho2 - minrho2) / (Nrho2 - 1)), maxrho2 + (30 * (maxrho2 - minrho2) / (Nrho2 - 1)), Nrho2 + 60);
    cosalpha_in2    = linspace(mincosalpha2 - (30 * (maxcosalpha2 - mincosalpha2) / (Ncosalpha2 - 1)), maxcosalpha2 + (30 * (maxcosalpha2 - mincosalpha2) / (Ncosalpha2 - 1)), Ncosalpha2 + 60);

    [RHO_IN2, COSALPHA_IN2] = meshgrid(rho_in2, cosalpha_in2);
    data_temp2 = data;
    data_temp2.RHO_OUT = RHO_IN2;
    data_temp2.COSALPHA_OUT = COSALPHA_IN2;

    data_temp2 = FFBP_POLAR_NFFT(data_temp2, center + 1, b, level + 1);

    data_temp2.Nearest  = data_temp2.Nearest    .* exp(-1j * (4 * pi / c * RHO_IN2) * fc);
    data_temp2.Linear   = data_temp2.Linear     .* exp(-1j * (4 * pi / c * RHO_IN2) * fc);
    data_temp2.Cubic    = data_temp2.Cubic      .* exp(-1j * (4 * pi / c * RHO_IN2) * fc);
    data_temp2.Spline   = data_temp2.Spline     .* exp(-1j * (4 * pi / c * RHO_IN2) * fc);
    data_temp2.knabNFFT = data_temp2.knabNFFT   .* exp(-1j * (4 * pi / c * RHO_IN2) * fc);
    data_temp2.rectNFFT = data_temp2.rectNFFT   .* exp(-1j * (4 * pi / c * RHO_IN2) * fc);
    data_temp2.knab     = data_temp2.knab       .* exp(-1j * (4 * pi / c * RHO_IN2) * fc);

    deltarho2       = mean(diff(RHO_IN2(1, :)));
    deltacosalpha2  = mean(diff(COSALPHA_IN2(:, 1)));

    % --- Interpolation of the image from the input local polar grid to the final global polar grid
    Reconstruction_Nearest2     = interp2(              RHO_IN2, COSALPHA_IN2, data_temp2.Nearest,  RHO_OUT2, COSALPHA_OUT2, 'nearest', ExtrapVal) .* exp(1j * (4 * pi / c * RHO_OUT2) * fc);
    Reconstruction_Linear2      = interp2(              RHO_IN2, COSALPHA_IN2, data_temp2.Linear,   RHO_OUT2, COSALPHA_OUT2, 'linear',  ExtrapVal) .* exp(1j * (4 * pi / c * RHO_OUT2) * fc);
    Reconstruction_Cubic2       = interp2(              RHO_IN2, COSALPHA_IN2, data_temp2.Cubic,    RHO_OUT2, COSALPHA_OUT2, 'cubic',   ExtrapVal) .* exp(1j * (4 * pi / c * RHO_OUT2) * fc);
    Reconstruction_Spline2      = interp2(              RHO_IN2, COSALPHA_IN2, data_temp2.Spline,   RHO_OUT2, COSALPHA_OUT2, 'spline',  ExtrapVal) .* exp(1j * (4 * pi / c * RHO_OUT2) * fc);
    Reconstruction_KnabNFFT2    = knabNFFT2D(ExtrapVal, RHO_IN2, COSALPHA_IN2, data_temp2.knabNFFT, RHO_OUT2, COSALPHA_OUT2, (1 / deltarho2) * 0.5, (1 / deltacosalpha2) * 0.5) .* exp(1j * (4 * pi / c * RHO_OUT2) * fc);
    Reconstruction_RectNFFT2    = rectNFFT2D(ExtrapVal, RHO_IN2, COSALPHA_IN2, data_temp2.rectNFFT, RHO_OUT2, COSALPHA_OUT2, (1 / deltarho2) * 0.5, (1 / deltacosalpha2) * 0.5) .* exp(1j * (4 * pi / c * RHO_OUT2) * fc);
    Reconstruction_Knab2        = knab2D(    ExtrapVal, RHO_IN2, COSALPHA_IN2, data_temp2.knab,     RHO_OUT2, COSALPHA_OUT2, (1 / deltarho2) * 0.5, (1 / deltacosalpha2) * 0.5) .* exp(1j * (4 * pi / c * RHO_OUT2) * fc);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % SUMMING UP THE RESULTS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    data.Nearest    = Reconstruction_Nearest1  + Reconstruction_Nearest2;
    data.Linear     = Reconstruction_Linear1   + Reconstruction_Linear2;
    data.Cubic      = Reconstruction_Cubic1    + Reconstruction_Cubic2;
    data.Spline     = Reconstruction_Spline1   + Reconstruction_Spline2;
    data.knabNFFT   = Reconstruction_KnabNFFT1 + Reconstruction_KnabNFFT2;
    data.rectNFFT   = Reconstruction_RectNFFT1 + Reconstruction_RectNFFT2;
    data.knab       = Reconstruction_Knab1     + Reconstruction_Knab2;

end

end
