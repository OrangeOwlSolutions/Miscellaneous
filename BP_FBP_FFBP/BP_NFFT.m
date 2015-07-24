function Reconstruction_BP = BP_NFFT(data, minInd, maxInd)

global c                                                        % --- Speed of light [m/s]
global fc                                                       % --- Center frequency [Hz]
global BW                                                       % --- Bandwidth [Hz]

% --- Image initialization with all zeros
Reconstruction_BP = zeros(1,numel(data.x_mat));

% --- Loop through every used pulse
for n = minInd : maxInd

    % --- Calculate the differential range for each pixel in the image [m]
    delta_r = sqrt((data.AntX(n)-data.x_mat).^2+(data.AntY(n)-data.y_mat).^2+(data.AntZ(n)-data.z_mat).^2)-data.R0(n);
    % --- Calculate the phase correction term
    PhaseCorrection = reshape(exp(1j*4*pi*fc*delta_r/c), 1, numel(delta_r));
    % --- Image update
    Reconstruction_BP = Reconstruction_BP + PhaseCorrection .* NFFT1_1D(data.phdata(:,n).', reshape(-2*BW/c*delta_r, 1, numel(delta_r)));

end

Reconstruction_BP = reshape(Reconstruction_BP, size(data.x_mat));

end

