function data = FFBP(data, a, b)

global c                                                        % --- Speed of light [m/s]

global fc                                                       % --- Center frequency [Hz]
global BW                                                       % --- Bandwidth [Hz]
global f1                                                       % --- Minimum frequency [Hz]
global f2                                                       % --- Maximum frequency [Hz]

center = a;
Dx = (data.AntX(center + 1) - data.AntX(center));
Dy = (data.AntY(center + 1) - data.AntY(center));
Dz = (data.AntZ(center + 1) - data.AntZ(center));
dx = Dx / sqrt(Dx^2 + Dy^2 + Dz^2);
dy = Dy / sqrt(Dx^2 + Dy^2 + Dz^2);
dz = Dz / sqrt(Dx^2 + Dy^2 + Dz^2);

% --- Creates the polar grid
RHO_OUT             = sqrt((data.x_mat - data.AntX(center)).^2 + (data.y_mat - data.AntY(center)).^2 + (data.z_mat - data.AntZ(center)).^2);
COSALPHA_OUT        = ((data.x_mat - data.AntX(center)) ./ RHO_OUT * dx + (data.y_mat - data.AntY(center)) ./ RHO_OUT * dy + (data.z_mat - data.AntZ(center)) ./ RHO_OUT * dz);
data.RHO_OUT        = RHO_OUT;
data.COSALPHA_OUT   = COSALPHA_OUT;

disp('Starting iteration number');
data.a = a;
data = FFBP_POLAR_NFFT(data, a, b, 0);
disp('Ending iteration number');

end
