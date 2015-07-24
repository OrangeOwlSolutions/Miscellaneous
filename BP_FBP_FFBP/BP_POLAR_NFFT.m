function data = BP_POLAR_NFFT(data, minInd, maxInd)

global c                                                        % --- Speed of light [m/s]
global fc                                                       % --- Center frequency [Hz]
global BW                                                       % --- Bandwidth [Hz]

% --- Partial image initialization with all zeros
data.im_final = zeros(1,numel(data.pho_mat));

% --- Loop through every used pulse
for n = minInd : maxInd

    % --- Calculate the differential range for each pixel in the image [m]
    xsi     = sqrt((data.AntX(n)-data.AntX(minInd)).^2+(data.AntY(n)-data.AntY(minInd)).^2+(data.AntZ(n)-data.AntZ(minInd)).^2);
    delta_r = sqrt(data.pho_mat.^2+xsi.^2-2*data.pho_mat.*xsi.*data.COSALPHA_OUT) - data.R0(n);
    % --- Calculate the phase correction term
    PhaseCorrection = reshape(exp(1j*4*pi*fc*delta_r/c), 1, numel(delta_r));
    % --- Image update
    data.im_final = data.im_final + PhaseCorrection.* NFFT1_1D(data.phdata(:, n).', reshape(-2 * BW / c * delta_r, 1, numel(delta_r)));

end

data.im_final = reshape(data.im_final,size(data.pho_mat));

end

