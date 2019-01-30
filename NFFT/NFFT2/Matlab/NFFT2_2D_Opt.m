function transf = NFFT2_2D_Opt(data, x1, x2, M1, M2, c, K, Max_Num_PFs, SBP_factor)
   
N=length(data);                                 % Number of input elements                              

% --- Algorithm parameters
alfa=(2-1/c)*pi-0.01;                           % Half-size of the support of the interpolation window
alfa_prime = ((2 - 1 / c) * pi);

mm = -K : K;
SBP_product = (SBP_factor * (2 * pi - pi / c) * K);

filename = strcat('Result_c_', num2str(c), '_NumPFs_', num2str(Max_Num_PFs), '_K_',num2str(K), '_SBPfactor_', num2str(SBP_factor), '.mat');
load(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPATIAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu1 = round(c * x1);                                          
kk1 = -M1 / 2 : M1 / 2 - 1;
xi1 = 2 * pi * kk1 / (c * M1);
[PSWFs, P, V_even , ~, K_Leg] = S0n(SBP_product, 2 * Max_Num_PFs, xi1 / (2 * pi - pi / c), 1e-30);
PSWFs = PSWFs.';
phi_cap1 = zeros(length(x1), length(mm));
[XX1, MM1] = meshgrid(x1, mm);
for p = 0 : 2 : K_Leg - 1
    temp = sqrt(pi ./ (2 .* alfa_prime * abs(c * XX1 - (MM1 + round(c * XX1))))) .* besselj(p + 0.5, alfa_prime * abs(c * XX1 - (MM1 + round(c * XX1))));
    indices = abs(c * XX1 - (MM1 + round(c * XX1))) == 0;
    temp(indices) = sqrt(pi) * (0.5 * alfa_prime * abs(c * XX1(indices) - (MM1(indices) + round(c * XX1(indices))))) .^ p ./ (2 * gamma(p + 1.5));
    phi_cap1 = phi_cap1 + alfa_prime * (sqrt(p + 0.5) / (2 * pi)) * sum(V_even(p / 2 + 1, 1 : length(c_opt)) .* c_opt) * ((2 * (-1i)^p .* temp)).';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPECTRAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Recovering the "optimal" spatial window by the Legendre Polynomials
phi1 = zeros(1, length(xi1));
for p = 0 : 2 : K_Leg - 1
    phi1 = phi1 + sum(V_even(p / 2 + 1, 1 : length(c_opt)) .* c_opt) * P(p + 1, :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPATIAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu2 = round(c * x2);                                         
kk2 = -M2 / 2 : M2 / 2 - 1;
xi2 = 2 * pi * kk2 / (c * M2);
[PSWFs, P, V_even , ~, K_Leg] = S0n(SBP_product, 2 * Max_Num_PFs, xi1 / (2 * pi - pi / c), 1e-30);
PSWFs = PSWFs.';
phi_cap2 = zeros(length(x2), length(mm));
[XX2, MM2] = meshgrid(x2, mm);
for p = 0 : 2 : K_Leg - 1
    temp = sqrt(pi ./ (2 .* alfa_prime * abs(c * XX2 - (MM2 + round(c * XX2))))) .* besselj(p + 0.5, alfa_prime * abs(c * XX2 - (MM2 + round(c * XX2))));
    indices = abs(c * XX2 - (MM2 + round(c * XX2))) == 0;
    temp(indices) = sqrt(pi) * (0.5 * alfa_prime * abs(c * XX2(indices) - (MM2(indices) + round(c * XX2(indices))))) .^ p ./ (2 * gamma(p + 1.5));
    phi_cap2 = phi_cap2 + alfa_prime * (sqrt(p + 0.5) / (2 * pi)) * sum(V_even(p / 2 + 1, 1 : length(c_opt)) .* c_opt) * ((2 * (-1i)^p .* temp)).';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPECTRAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Recovering the "optimal" spatial window by the Legendre Polynomials
phi2 = zeros(1, length(xi2));
for p = 0 : 2 : K_Leg - 1
    phi2 = phi2 + sum(V_even(p/2+1,1:length(c_opt)).*c_opt) * P(p+1,:);
end

[PHI2, PHI1] = meshgrid(phi2, phi1);

%%%%%%%%%%
% STEP 1 %
%%%%%%%%%%
u = zeros(c * M1, c * M2);
for k = 1 : N
    for l1 = -K : K
        for l2 = -K : K,
            p1 = mod(mu1(k) + l1 + M1 * c / 2, c * M1) + 1;
            p2 = mod(mu2(k) + l2 + M2 * c / 2, c * M2) + 1;
            u(p1, p2) = u(p1, p2) + data(k) * phi_cap1(k, l1 + K + 1) * phi_cap2(k, l2 + K + 1);
        end
    end
end

%%%%%%%%%%
% STEP 2 %
%%%%%%%%%%
% FFT_Mc=M1*M2*(c*c)*ifftshift(ifft2(fftshift(u)));
FFT_Mc = fftshift(fft2(ifftshift(u)));
% FFT_Mc=fft2(u);

%%%%%%%%%%
% STEP 3 %
%%%%%%%%%%
transf = FFT_Mc(1 + (c - 1) * M1 / 2 : (c + 1) * M1 / 2, 1 + (c - 1) * M2 / 2 : (c + 1) * M2 / 2) ./ (PHI1 .* PHI2);

end 
