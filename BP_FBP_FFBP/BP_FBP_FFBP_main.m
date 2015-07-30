clear all; 
close all;
clc;

global c fc f1 f2 BW NumFFTSamples maxWr

c          = 299792458;                                          % --- Speed of light [m/s]

%%%%%%%%%%%%%%%%%%
% SAR PARAMETERS %
%%%%%%%%%%%%%%%%%%
data.BW         = 600*1e6;                                       % --- Bandwidth [Hz]
BW              = 600*1e6;                                       % --- Bandwidth [Hz]
fc              = 10*1e9;                                        % --- Center frequency [Hz]
%% --- ??
data.K          = 216;                                           % --- Number of frequency samples per pulses
%% --- ??
data.minF       = (fc - data.BW/2.0);                            % --- Minimum frequency for the illuminating pulse [Hz]
%% --- ??
data.deltaF     = data.BW /data.K;                               % --- Frequency step size [Hz]

%% --- ??
Np              = 128;                                           % --- Number of pulses for the antenna aperture
% Np              = 10;                                           % --- Number of pulses for the antenna aperture
R               = 10*1e3;                                        % --- Range distance to the scene center [m]
depAngle        = pi/6;                                          % --- Depression angle [m]
maxAz           = 3*pi/180;                                      % --- Maximum angular size of the aperture [radians]

Wx              = 10.24;                                         % --- Scene dimension along the x axis 
Wy              = 10.24;                                         % --- Scene dimension along the y axis

NumFFTSamples   = 4096;                                          % --- Length of the FFT 
data.Nfft       = 4096;                                          % --- Length of the FFT 

f1 = data.minF(1);                                               % --- Minimum frequency
f2 = data.minF(1) + data.K *data.deltaF;                         % --- Maximum frequency

%%%%%%%%%%%%%%%%
% SAR GEOMETRY %
%%%%%%%%%%%%%%%%
AntAzimuth  = (0:Np-1)/Np*maxAz-maxAz/2;                         % --- Angular position of the antenna for every pulses

data.AntX   = R*cos(depAngle)*ones(1,Np);                        % --- x-coordinates of the antenna for every pulse
data.AntY   = R*cos(depAngle)*tan(AntAzimuth);                   % --- y-coordinates of the antenna for every pulse
data.AntZ   = R*sin(depAngle)*ones(1,Np);                        % --- z-coordinates of the antenna for every pulse
fprintf('Aperture length: %f m\n', max(data.AntY)-min(data.AntY));

AntAzimuth  = unwrap(atan2(data.AntY,data.AntX));

maxAz       = max(AntAzimuth) - min(AntAzimuth);                 % --- Maximum angular size of the aperture [radians]

data.deltaAz= abs(mean(diff(AntAzimuth)));                       % --- Average azimuth angle step size (radians)

% --- Distances from antenna position for every pulse and scene center
data.R0 = sqrt((data.AntX).^2+(data.AntY).^2+(data.AntZ).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAR PERFORMANCE PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Maximum scene size
maxWr           = c/(2*data.deltaF);                            % --- Maximum scene size of the image in range [m]
data.maxWr      = c/(2*data.deltaF);                            % --- Maximum scene size of the image in range [m]
data.maxWx      = c/(2*data.deltaAz*data.minF);                 % --- Maximum scene size of the image in azimuth [m] 
fprintf('Maximum Scene Size: %.2f m in range, %.2f m in cross-range\n', data.maxWr, data.maxWx);

% --- Resolutions
data.dr         = c/(2*data.deltaF*data.K);                     % --- Resolution in range [m]
data.dx         = c/(2*maxAz*data.minF);                        % --- Resolution in azimuth [m]
fprintf('Resolution: %.2fm in range, %.2f m in cross-range\n', data.dr, data.dx);

%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTING UP IMAGE GRID %
%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- ??
Xstep = 0.08;                                                   % --- Sampling step of the image along x [m]
Ystep = 0.08;                                                   % --- Sampling step of the image along y [m]

Nx    = Wx/Xstep;                                               % --- Number of image points along x
Ny    = Wy/Ystep;                                               % --- Number of image points along y

x_vec = linspace(-Wx/2,Wx/2,Nx);
y_vec = linspace(-Wy/2,Wy/2,Ny);

[data.x_mat,data.y_mat] = meshgrid(x_vec,y_vec);
data.z_mat = 0*ones(size(data.x_mat));

%%%%%%%%%%%%%%%%%%%
% DATA GENERATION %
%%%%%%%%%%%%%%%%%%%
x = 0;
y = 0;
z = 0;

data.phdata = ones(data.K,Np);
% --- Loop over the pulses
for n=1:Np
    % --- Loop over the frequencies
    for k=1:data.K
        fk = data.minF+(k-1)*data.deltaF;
        deltaR = sqrt((data.AntX(n)-x).^2+(data.AntY(n)-y).^2+(data.AntZ(n)-z).^2)-data.R0(n);
        data.phdata(k,n) = 1*exp((-1j*4*pi*fk*deltaR)/c);       
    end
end

%%%%%%%%%%%%%%%%%%%
% PLOT PARAMETERS %
%%%%%%%%%%%%%%%%%%%
dyn_range = 40;                     % --- Dynamic range for image plotting
FontSize  = 14;                     % --- Font size

a = 1; %first used pulse 
b = 128; %last used pulse

% a = 1; %first used pulse 
% b = 4; %last used pulse

%%%%%%%%%%%%%%%%%
% BP USING NFFT %
%%%%%%%%%%%%%%%%%
% Reconstruction_BP = BP_NFFT(data, a, b);
% DisplayImage(1, x_vec, y_vec, Reconstruction_BP, dyn_range, 'BP - NFFT', FontSize)

%%%%%%%%%%%%%%%%%%%%%%%%%%
% BP USING INTERPOLATORS %
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction_BP_Interp = BP(data,a,b);
% DisplayImage(2, x_vec, y_vec, Reconstruction_BP_Interp.Nearest,  dyn_range, 'BP - Nearest',                 FontSize)
% DisplayImage(3, x_vec, y_vec, Reconstruction_BP_Interp.Linear,   dyn_range, 'BP - Linear',                  FontSize)
% DisplayImage(4, x_vec, y_vec, Reconstruction_BP_Interp.Cubic,    dyn_range, 'BP - Cubic',                   FontSize)
% DisplayImage(5, x_vec, y_vec, Reconstruction_BP_Interp.Spline,   dyn_range, 'BP - Spline',                  FontSize)
% DisplayImage(6, x_vec, y_vec, Reconstruction_BP_Interp.Knab,     dyn_range, 'BP - Knab window',             FontSize)
% DisplayImage(7, x_vec, y_vec, Reconstruction_BP_Interp.RectNFFT, dyn_range, 'BP - NFFT with rect window',   FontSize)
% DisplayImage(8, x_vec, y_vec, Reconstruction_BP_Interp.KnabNFFT, dyn_range, 'BP - NFFT with Knab window',   FontSize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FBP USING INTERPOLATORS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction_FBP_Interp = FBP(data, a, b);
% DisplayImage(9,  x_vec, y_vec, Reconstruction_FBP_Interp.Nearest,   dyn_range, 'FBP - Nearest',                 FontSize)
% DisplayImage(10, x_vec, y_vec, Reconstruction_FBP_Interp.Linear,    dyn_range, 'FBP - Linear',                  FontSize)
% DisplayImage(11, x_vec, y_vec, Reconstruction_FBP_Interp.Cubic,     dyn_range, 'FBP - Cubic',                   FontSize)
% DisplayImage(12, x_vec, y_vec, Reconstruction_FBP_Interp.Spline,    dyn_range, 'FBP - Spline',                  FontSize)
% DisplayImage(13, x_vec, y_vec, Reconstruction_FBP_Interp.knab,      dyn_range, 'FBP - Knab window',             FontSize)
% DisplayImage(14, x_vec, y_vec, Reconstruction_FBP_Interp.knabNFFT,  dyn_range, 'FBP - NFFT with Knab window',   FontSize)
% DisplayImage(15, x_vec, y_vec, Reconstruction_FBP_Interp.rectNFFT,  dyn_range, 'FBP - NFFT with rect window',   FontSize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFBP USING INTERPOLATORS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction_FFBP_Interp = FFBP(data, a, b);
% DisplayImage(16, x_vec, y_vec, Reconstruction_FFBP_Interp.Nearest,     dyn_range, 'FFBP - Nearest',                FontSize)
% DisplayImage(17, x_vec, y_vec, Reconstruction_FFBP_Interp.Linear,      dyn_range, 'FFBP - Linear',                 FontSize)
% DisplayImage(18, x_vec, y_vec, Reconstruction_FFBP_Interp.Cubic,       dyn_range, 'FFBP - Cubic',                  FontSize)
% DisplayImage(19, x_vec, y_vec, Reconstruction_FFBP_Interp.Spline,      dyn_range, 'FFBP - Spline',                 FontSize)
% DisplayImage(20, x_vec, y_vec, Reconstruction_FFBP_Interp.knab,        dyn_range, 'FFBP - Knab window',            FontSize)
% DisplayImage(21, x_vec, y_vec, Reconstruction_FFBP_Interp.knabNFFT,    dyn_range, 'FFBP - NFFT with Knab window',  FontSize)
% DisplayImage(22, x_vec, y_vec, Reconstruction_FFBP_Interp.rectNFFT,    dyn_range, 'FFBP - NFFT with rect window',  FontSize)

load Test.mat
load Test2.mat

% return

%%%%%%%%%%%%%%%%%%
% NORMALIZATIONS %
%%%%%%%%%%%%%%%%%%
Reference_reconstruction = (Reconstruction_BP) ./ max(max(abs(Reconstruction_BP)));

Reconstruction_BP_Nearest       = Reconstruction_BP_Interp.Nearest      ./ max(max(abs(Reconstruction_BP_Interp.Nearest)));
Reconstruction_BP_Linear        = Reconstruction_BP_Interp.Linear       ./ max(max(abs(Reconstruction_BP_Interp.Linear)));
Reconstruction_BP_Cubic         = Reconstruction_BP_Interp.Cubic        ./ max(max(abs(Reconstruction_BP_Interp.Cubic)));
Reconstruction_BP_Spline        = Reconstruction_BP_Interp.Spline       ./ max(max(abs(Reconstruction_BP_Interp.Spline)));
Reconstruction_BP_Knab          = Reconstruction_BP_Interp.Knab         ./ max(max(abs(Reconstruction_BP_Interp.Knab)));
Reconstruction_BP_RectNFFT      = Reconstruction_BP_Interp.RectNFFT     ./ max(max(abs(Reconstruction_BP_Interp.RectNFFT)));
Reconstruction_BP_KnabNFFT      = Reconstruction_BP_Interp.KnabNFFT     ./ max(max(abs(Reconstruction_BP_Interp.KnabNFFT)));

Reconstruction_FBP_Nearest      = Reconstruction_FBP_Interp.Nearest     ./ max(max(abs(Reconstruction_FBP_Interp.Nearest)));
Reconstruction_FBP_Linear       = Reconstruction_FBP_Interp.Linear      ./ max(max(abs(Reconstruction_FBP_Interp.Linear)));
Reconstruction_FBP_Cubic        = Reconstruction_FBP_Interp.Cubic       ./ max(max(abs(Reconstruction_FBP_Interp.Cubic)));
Reconstruction_FBP_Spline       = Reconstruction_FBP_Interp.Spline      ./ max(max(abs(Reconstruction_FBP_Interp.Spline)));
Reconstruction_FBP_Knab         = Reconstruction_FBP_Interp.knab        ./ max(max(abs(Reconstruction_FBP_Interp.knab)));
Reconstruction_FBP_RectNFFT     = Reconstruction_FBP_Interp.rectNFFT    ./ max(max(abs(Reconstruction_FBP_Interp.rectNFFT)));
Reconstruction_FBP_KnabNFFT     = Reconstruction_FBP_Interp.knabNFFT    ./ max(max(abs(Reconstruction_FBP_Interp.knabNFFT)));

Reconstruction_FFBP_Nearest     = Reconstruction_FFBP_Interp.Nearest    ./ max(max(abs(Reconstruction_FFBP_Interp.Nearest)));
Reconstruction_FFBP_Linear      = Reconstruction_FFBP_Interp.Linear     ./ max(max(abs(Reconstruction_FFBP_Interp.Linear)));
Reconstruction_FFBP_Cubic       = Reconstruction_FFBP_Interp.Cubic      ./ max(max(abs(Reconstruction_FFBP_Interp.Cubic)));
Reconstruction_FFBP_Spline      = Reconstruction_FFBP_Interp.Spline     ./ max(max(abs(Reconstruction_FFBP_Interp.Spline)));
Reconstruction_FFBP_Knab        = Reconstruction_FFBP_Interp.knab       ./ max(max(abs(Reconstruction_FFBP_Interp.knab)));
Reconstruction_FFBP_RectNFFT    = Reconstruction_FFBP_Interp.rectNFFT   ./ max(max(abs(Reconstruction_FFBP_Interp.rectNFFT)));
Reconstruction_FFBP_KnabNFFT    = Reconstruction_FFBP_Interp.knabNFFT   ./ max(max(abs(Reconstruction_FFBP_Interp.knabNFFT)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOT MEAN SQUARE ERRORS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Reconstruction_BP_Nearest       = Reconstruction_BP_Interp.Nearest      ./ max(max(abs(Reconstruction_BP_Interp.Nearest)));
Reconstruction_BP_Linear        = Reconstruction_BP_Interp.Linear       ./ max(max(abs(Reconstruction_BP_Interp.Linear)));
Reconstruction_BP_Cubic         = Reconstruction_BP_Interp.Cubic        ./ max(max(abs(Reconstruction_BP_Interp.Cubic)));
Reconstruction_BP_Spline        = Reconstruction_BP_Interp.Spline       ./ max(max(abs(Reconstruction_BP_Interp.Spline)));
Reconstruction_BP_Knab          = Reconstruction_BP_Interp.Knab         ./ max(max(abs(Reconstruction_BP_Interp.Knab)));
Reconstruction_BP_RectNFFT      = Reconstruction_BP_Interp.RectNFFT     ./ max(max(abs(Reconstruction_BP_Interp.RectNFFT)));
Reconstruction_BP_KnabNFFT      = Reconstruction_BP_Interp.KnabNFFT     ./ max(max(abs(Reconstruction_BP_Interp.KnabNFFT)));

RMS.BP_Nearest      = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_BP_Nearest )).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.BP_Linear       = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_BP_Linear  )).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.BP_Cubic        = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_BP_Cubic   )).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.BP_Spline       = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_BP_Spline  )).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.BP_Knab         = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_BP_Knab    )).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.BP_RectNFFT     = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_BP_RectNFFT)).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.BP_KnabNFFT     = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_BP_KnabNFFT)).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));

RMS.FBP_Nearest      = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_FBP_Nearest )).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.FBP_Linear       = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_FBP_Linear  )).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.FBP_Cubic        = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_FBP_Cubic   )).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.FBP_Spline       = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_FBP_Spline  )).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.FBP_Knab         = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_FBP_Knab    )).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.FBP_RectNFFT     = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_FBP_RectNFFT)).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.FBP_KnabNFFT     = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_FBP_KnabNFFT)).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));

RMS.FFBP_Nearest      = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_FFBP_Nearest )).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.FFBP_Linear       = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_FFBP_Linear  )).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.FFBP_Cubic        = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_FFBP_Cubic   )).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.FFBP_Spline       = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_FFBP_Spline  )).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.FFBP_Knab         = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_FFBP_Knab    )).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.FFBP_RectNFFT     = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_FFBP_RectNFFT)).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));
RMS.FFBP_KnabNFFT     = 100 * sqrt((sum(sum(abs(((Reference_reconstruction - Reconstruction_FFBP_KnabNFFT)).^2)))) / sum(sum(abs((Reference_reconstruction).^2))));

%%%%%%%%%%%%%%%%%
% ERROR FIGURES %
%%%%%%%%%%%%%%%%%

LineWidth = 2;

figure
plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_BP_Nearest(:, Np/2))),  'g', 'LineWidth', LineWidth)
hold on
plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_BP_Linear(:, Np/2))),   'c', 'LineWidth', LineWidth)
plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_BP_Cubic(:, Np/2))),    'r', 'LineWidth', LineWidth)
plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_BP_Spline(:, Np/2))),   'k', 'LineWidth', LineWidth)
% plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_BP_RectNFFT(:, Np/2))), 'm', 'LineWidth', LineWidth)
plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_BP_KnabNFFT(:, Np/2))), 'b', 'LineWidth', LineWidth)
% plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_BP_Knab(:, Np/2))),     'y', 'LineWidth', LineWidth)
h=title('BP');
set(h, 'FontSize', FontSize);
h = xlabel('x [m]');
set(h, 'FontSize', FontSize);
h = ylabel('Error [dB]');
set(h, 'FontSize', FontSize);
set(gca, 'FontSize', FontSize);
% legend('Nearest','Linear','Cubic','Spline','Rect NUFFT','Knab NUFFT','Knab')
legend('Nearest','Linear','Cubic','Spline','Knab NUFFT')
axis([min(x_vec) max(x_vec) -300 0])

figure
plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_FBP_Nearest(:, Np/2))),  'g', 'LineWidth', LineWidth)
hold on
plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_FBP_Linear(:, Np/2))),   'c', 'LineWidth', LineWidth)
plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_FBP_Cubic(:, Np/2))),    'r', 'LineWidth', LineWidth)
plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_FBP_Spline(:, Np/2))),   'k', 'LineWidth', LineWidth)
% plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_FBP_RectNFFT(:, Np/2))), 'm', 'LineWidth', LineWidth)
plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_FBP_KnabNFFT(:, Np/2))), 'b', 'LineWidth', LineWidth)
% plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_FBP_Knab(:, Np/2))),     'y', 'LineWidth', LineWidth)
h=title('FBP');
set(h, 'FontSize', FontSize);
h = xlabel('x [m]');
set(h, 'FontSize', FontSize);
h = ylabel('Error [dB]');
set(h, 'FontSize', FontSize);
set(gca, 'FontSize', FontSize);
% legend('Nearest','Linear','Cubic','Spline','Rect NUFFT','Knab NUFFT','Knab')
legend('Nearest','Linear','Cubic','Spline','Knab NUFFT')
axis([min(x_vec) max(x_vec) -250 0])

figure
plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_FFBP_Nearest(:, Np/2))),  'g', 'LineWidth', LineWidth)
hold on
plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_FFBP_Linear(:, Np/2))),   'c', 'LineWidth', LineWidth)
plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_FFBP_Cubic(:, Np/2))),    'r', 'LineWidth', LineWidth)
plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_FFBP_Spline(:, Np/2))),   'k', 'LineWidth', LineWidth)
% plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_FFBP_RectNFFT(:, Np/2))), 'm', 'LineWidth', LineWidth)
plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_FFBP_KnabNFFT(:, Np/2))), 'b', 'LineWidth', LineWidth)
% plot(x_vec, 20*log10(abs(Reference_reconstruction(:, Np/2) - Reconstruction_FFBP_Knab(:, Np/2))),     'y', 'LineWidth', LineWidth)
h=title('FBP');
set(h, 'FontSize', FontSize);
h = xlabel('x [m]');
set(h, 'FontSize', FontSize);
h = ylabel('Error [dB]');
set(h, 'FontSize', FontSize);
set(gca, 'FontSize', FontSize);
% legend('Nearest','Linear','Cubic','Spline','Rect NUFFT','Knab NUFFT','Knab')
legend('Nearest','Linear','Cubic','Spline','Knab NUFFT')
axis([min(x_vec) max(x_vec) -150 0])

