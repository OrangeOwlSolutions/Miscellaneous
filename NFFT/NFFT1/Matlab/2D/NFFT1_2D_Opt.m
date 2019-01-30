function result = NFFT1_2D_Opt(data, x1, x2, c, K, Max_Num_PFs, SBP_factor)

N1 = size(data, 1);
N2 = size(data, 2);

% --- Algorithm parameters
alfa_prime = ((2 - 1 / c) * pi);

filename = strcat('Result_c_', num2str(c), '_NumPFs_', num2str(Max_Num_PFs), '_K_',num2str(K), '_SBPfactor_', num2str(SBP_factor), '.mat');
load(filename)

mm = -K : K;
SBP_product = (SBP_factor * (2 * pi - pi / c) * K);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPATIAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk1 = -N1 / 2 : N1 / 2 - 1;
xi1 = 2 * pi * kk1 / (c * N1);

[PSWFs, P, V_even , ~, K_Leg] = S0n(SBP_product, 2 * Max_Num_PFs, xi1 / (2 * pi - pi / c), 1e-30);
PSWFs = PSWFs.';

% --- Recovering the "optimal" spatial window by the Legendre Polynomials
phi1 = zeros(1, length(xi1));
for p = 0 : 2 : K_Leg - 1
    phi1 = phi1 + sum(V_even(p / 2 + 1, 1 : length(c_opt)) .* c_opt) * P(p + 1, :);
end

kk2 = -N2 / 2 : N2 / 2 - 1;
xi2 = 2 * pi * kk2 / (c * N2);

[PSWFs, P, V_even , ~, K_Leg] = S0n(SBP_product, 2 * Max_Num_PFs, xi2 / (2 * pi - pi / c), 1e-30);
PSWFs = PSWFs.';

% --- Recovering the "optimal" spatial window by the Legendre Polynomials
phi2 = zeros(1, length(xi2));
for p = 0 : 2 : K_Leg - 1
    phi2 = phi2 + sum(V_even(p/2+1,1:length(c_opt)).*c_opt) * P(p+1,:);
end

[PHI2, PHI1] = meshgrid(phi2, phi1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPECTRAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu1 = round(c * x1);  

phi_cap1 = zeros(length(x1), length(mm));
[XX1, MM1] = meshgrid(x1, mm);
for p = 0 : 2 : K_Leg - 1
    temp = sqrt(pi ./ (2 .* alfa_prime * abs(c * XX1 - (MM1 + round(c * XX1))))) .* besselj(p + 0.5, alfa_prime * abs(c * XX1 - (MM1 + round(c * XX1))));
    indices = abs(c * XX1 - (MM1 + round(c * XX1))) == 0;
    temp(indices) = sqrt(pi) * (0.5 * alfa_prime * abs(c * XX1(indices) - (MM1(indices) + round(c * XX1(indices))))) .^ p ./ (2 * gamma(p + 1.5));
    phi_cap1 = phi_cap1 + alfa_prime * (sqrt(p + 0.5) / (2 * pi)) * sum(V_even(p / 2 + 1, 1 : length(c_opt)) .* c_opt) * ((2 * (-1i)^p .* temp)).';
end

mu2 = round(c * x2);  

phi_cap2 = zeros(length(x2), length(mm));
[XX2, MM2] = meshgrid(x2, mm);
for p = 0 : 2 : K_Leg - 1
    temp = sqrt(pi ./ (2 .* alfa_prime * abs(c * XX2 - (MM1 + round(c * XX2))))) .* besselj(p + 0.5, alfa_prime * abs(c * XX2 - (MM2 + round(c * XX2))));
    indices = abs(c * XX2 - (MM2 + round(c * XX2))) == 0;
    temp(indices) = sqrt(pi) * (0.5 * alfa_prime * abs(c * XX2(indices) - (MM2(indices) + round(c * XX2(indices))))) .^ p ./ (2 * gamma(p + 1.5));
    phi_cap2 = phi_cap2 + alfa_prime * (sqrt(p + 0.5) / (2 * pi)) * sum(V_even(p / 2 + 1, 1 : length(c_opt)) .* c_opt) * ((2 * (-1i)^p .* temp)).';
end

%%%%%%%%%%
% STEP 1 %
%%%%%%%%%%
u = zeros(c * N1, c * N2);
u((c - 1) * N1 / 2 + 1 : (c + 1) * N1 / 2, (c - 1) * N2 / 2 + 1 : (c + 1) * N2 / 2) = data ./ (PHI1 .* PHI2);

%%%%%%%%%%
% STEP 2 %
%%%%%%%%%%
U = fft2(ifftshift(u));

%%%%%%%%%%
% STEP 3 %
%%%%%%%%%%
N = length(x1);
result = zeros(1, N);
for l = 1 : N,
    for m1 = -K : K,
        result(l) = result(l) + sum(phi_cap1(l, m1 + K + 1) * phi_cap2(l, 1 : 2 * K + 1).*...
                                                      U(mod(mu1(l) + m1, c * N1) + 1, mod(mu2(l) + (-K : K), c * N2) + 1));
    end
end

end
