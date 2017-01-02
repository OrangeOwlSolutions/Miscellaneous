function F = NUFFT_3T_2D_LG_v3(x, y, f, s, t, eps)

Maxx = max(x); minx = min(x); xb = (Maxx + minx) / 2; X1 = abs(Maxx - xb);
Maxy = max(y); miny = min(y); yb = (Maxy + miny) / 2; Y1 = abs(Maxy - yb);
Maxs = max(s); mins = min(s); sb = (Maxs + mins) / 2; S  = abs(Maxs - sb);
Maxt = max(t); mint = min(t); tb = (Maxt + mint) / 2; T  = abs(Maxt - tb);

% --- Precision dependent parameters
R = 2.2;
msp = (2 + 2 * R^2 * (-log(eps / 76)) / (pi * (R^2 - 2)));

% --- Dx Ds
Mrx  = 2 * ceil((X1 * S / pi) * R^2 + R * msp);    
kk   = X1 / (pi / 2);
x    = x / kk;  
s    = s * kk;
Maxx = Maxx / kk;
minx = minx / kk;
Maxs = Maxs * kk;
mins = mins * kk;
xb   = (Maxx + minx) / 2;
sb   = (Maxs + mins) / 2;
Dx   = 2 * pi / Mrx;
Ds   = 1;

% --- Dy Dt
Mry  = 2 * ceil((Y1 * T / pi) * R^2 + R * msp);
kk   = Y1 / (pi / 2);
y    = y / kk;  
t    = t * kk;
Maxy = Maxy / kk;
miny = miny / kk;
Maxt = Maxt * kk;
mint = mint * kk;
yb   = (Maxy + miny) / 2;
tb   = (Maxt + mint) / 2;
Dy   = 2 * pi / Mry;
Dt   = 1;

% --- Convolution constants
b   = (msp - 2) / (4 * pi);  
t1  = 1 / (4 * b);
msp = floor(msp - 2)
np  = -msp : msp;
E3  = exp(-t1 * np.^2);

% --- Coefficients
c   = exp(-1i * sb * x) .* exp(-1i * tb * y) .* f;

% --- 1) Convolution: f_tau
Nxy   = length(x);
f_tau = zeros(Mry, Mrx); 
temp  = f_tau;
for kk = 1 : Nxy,
    
    % --- x coordinate
    nxn   = floor(Mrx / 2 + (x(kk) - xb) / Dx);
    difx  = (Mrx / 2 + (x(kk) - xb) / Dx) - nxn;
    
    % --- y coordinate
    nyn   = floor(Mry / 2 + (y(kk) - yb) / Dy);
    dify  = (Mry / 2 + (y(kk) - yb) / Dy) - nyn;
    
    % --- Create Gaussian 2D submatrix
    temp1 = c(kk) * exp(-t1 * difx^2) * exp(2 * t1 * difx * np) .* E3;
    temp2 =         exp(-t1 * dify^2) * exp(2 * t1 * dify * np) .* E3;
    temp(nyn + 1 - msp : nyn + 1 + msp, nxn + 1 - msp : nxn + 1 + msp) = (temp2.') * temp1;
    
    f_tau = temp + f_tau;
    temp  = zeros(Mry, Mrx);  
end

% --- 2) Compensation: f_st
ns1 = -Mrx / 2 : Mrx / 2 - 1;
ns2 = -Mry / 2 : Mry / 2 - 1;
f_st = ((1 / (4 * pi * b))^2) * f_tau .* ((exp((b * Dy^2) * ns2.^2).') * exp((b * Dx^2) * ns1.^2)); 

% --- 3) Transform: F_st    
F_st = ifftshift(fft2(fftshift(f_st)));

% --- 4) Convolution: F_t
Ns  = length(s);
F_t = zeros(1, Ns);
for kk = 1 : Ns
    
    % --- s coordinate
    nsn  = floor(Mrx / 2 + (s(kk) - sb) / Ds);
    difs = (Mrx / 2 + (s(kk) - sb) / Ds) - nsn;
    
    % --- t coordinate
    ntn  = floor(Mry / 2 + (t(kk) - tb) / Dt);
    dift = (Mry / 2 + (t(kk) - tb) / Dt) - ntn;
    
    % --- Create Gaussian 2D submatrix
    temp3  = ((exp(-t1 * dift^2) * exp(2 * t1 * dift * np) .* E3).') * (exp(-t1 * difs^2) * exp(2 * t1 * difs * np) .* E3);
    F_t(kk)= sum(sum(F_st(ntn + 1 - msp : ntn + 1 + msp, nsn + 1 - msp : nsn + 1 + msp) .* temp3));      
end

% --- 5) Compensation: F
F = F_t .* exp(-1i * (s - sb) * xb) .* exp(b * Dx^2 * (s - sb).^2);
F = F   .* exp(-1i * (t - tb) * yb) .* exp(b * Dy^2 * (t - tb).^2);

end
