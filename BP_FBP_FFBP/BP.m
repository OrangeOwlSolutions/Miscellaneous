function Reconstruction = BP(data, minInd, maxInd)

global c                                                        % --- Speed of light [m/s]
global NumFFTSamples                                            % --- Length of the FFT
global maxWr                                                    % --- Maximum scene size of the image in range [m]

% --- Range to every bin in the range profile [m]
data.r_vec = linspace(-NumFFTSamples/2, NumFFTSamples/2-1, NumFFTSamples) * maxWr / NumFFTSamples;

% --- Reconstruction initialization with all zeros
Reconstruction.KnabNFFT     = zeros(size(data.x_mat));
Reconstruction.Spline       = zeros(size(data.x_mat));
Reconstruction.Linear       = zeros(size(data.x_mat));
Reconstruction.Cubic        = zeros(size(data.x_mat));
Reconstruction.Nearest      = zeros(size(data.x_mat));
Reconstruction.Knab         = zeros(size(data.x_mat));
Reconstruction.RectNFFT     = zeros(size(data.x_mat));

% --- Loop through every used pulse
for n = minInd : maxInd
    
    % --- Zero padded raw data
    rc = fftshift(ifft(data.phdata(:,n), NumFFTSamples));
    % --- Calculate the differential range for each pixel in the image [m]
    delta_r = sqrt((data.AntX(n)-data.x_mat).^2+(data.AntY(n)-data.y_mat).^2+(data.AntZ(n)-data.z_mat).^2)-data.R0(n);
    % --- Calculate the phase correction term
    PhaseCorrection = exp(1i*4*pi*data.minF/c*delta_r);
    % --- Determine which pixels fall within the range swath
    I = find(and(delta_r > min(data.r_vec), delta_r < max(data.r_vec)));

    % --- Image update
    Reconstruction.Nearest(I)  = Reconstruction.Nearest(I)  + interp1(data.r_vec,rc,  delta_r(I), 'nearest')                                                                                                     .* PhaseCorrection(I);
    Reconstruction.Linear(I)   = Reconstruction.Linear(I)   + interp1(data.r_vec,rc,  delta_r(I), 'linear')                                                                                                      .* PhaseCorrection(I);
    Reconstruction.Cubic(I)    = Reconstruction.Cubic(I)    + interp1(data.r_vec,rc,  delta_r(I), 'cubic')                                                                                                       .* PhaseCorrection(I);
    Reconstruction.Spline(I)   = Reconstruction.Spline(I)   + interp1(data.r_vec,rc,  delta_r(I), 'spline')                                                                                                      .* PhaseCorrection(I);

    Reconstruction.KnabNFFT(I) = Reconstruction.KnabNFFT(I) + knabNFFT1D(data.r_vec.', rc.', delta_r(I).', 0.5*NumFFTSamples/maxWr).'                                                                              .* PhaseCorrection(I);
    Reconstruction.RectNFFT(I) = Reconstruction.RectNFFT(I) + rectNFFT1D(data.r_vec.', rc.', delta_r(I).', 0.5*NumFFTSamples/maxWr).'                                                                              .* PhaseCorrection(I);
    Reconstruction.Knab(I)     = Reconstruction.Knab(I)     + reshape(knab1D(data.r_vec,rc.', delta_r(I),                0.5*NumFFTSamples/maxWr), size(Reconstruction.Knab(I),1), size(Reconstruction.Knab(I),2)) .* PhaseCorrection(I);

end

return
