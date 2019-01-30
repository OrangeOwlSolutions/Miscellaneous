function result = NFFT1_1D_Opt(data, x, c, K, Max_Num_PFs, SBP_factor)

N = length(data);

% --- Algorithm parameters
alfa_prime = ((2 - 1 / c) * pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPATIAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk = -N / 2 : N / 2 - 1;
xi = 2 * pi * kk / (c * N);

%%%%%%%
% PFs %
%%%%%%%
SBP_product = (SBP_factor * (2 * pi - pi / c) * K);
[PSWFs, P, V_even1, ~, K_Leg] = S0n(SBP_product, 2 * Max_Num_PFs, xi / (2 * pi - pi / c), 1e-30);
PSWFs = PSWFs.';

filename = strcat('Result_c_', num2str(c), '_NumPFs_', num2str(Max_Num_PFs), '_K_',num2str(K), '_SBPfactor_', num2str(SBP_factor), '.mat');
load(filename)

% --- Recovering the "optimal" spatial window 
phi = zeros(1, length(xi));
for p = 0 : 2 : K_Leg - 1
    phi = phi + sum(V_even(p / 2 + 1, 1 : length(c_opt)) .* c_opt) * P(p + 1, :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPECTRAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = round(c * x);  
mm = -K : K;

[KK, MU] = meshgrid(mm, mu);

spectrum_phi = zeros(length(x), length(mm));
[XX, MM] = meshgrid(x, mm);
for p=0 : 2 : K_Leg - 1
    temp = sqrt(pi ./ (2 .* alfa_prime * abs(c * XX - (MM + round(c * XX))))) .* besselj(p + 0.5, alfa_prime * abs(c * XX - (MM + round(c * XX))));
    indices = abs(c * XX - (MM + round(c * XX))) == 0;
    temp(indices) = sqrt(pi) * (0.5 * alfa_prime * abs(c * XX(indices) - (MM(indices) + round(c * XX(indices))))) .^ p ./ (2 * gamma(p + 1.5));
    spectrum_phi = spectrum_phi + alfa_prime * (sqrt(p + 0.5) / (2 * pi)) * sum(V_even(p / 2 + 1, 1 : length(c_opt)) .* c_opt) * ((2 * (-1i)^p .* temp)).';
end

%%%%%%%%%%
% STEP 1 %
%%%%%%%%%%
u = zeros(1, c * N);
u((c - 1) * N / 2 + 1 : (c + 1) * N /2) = data ./ phi;

%%%%%%%%%%
% STEP 2 %
%%%%%%%%%%
U=fft(ifftshift(u));

%%%%%%%%%%
% STEP 3 %
%%%%%%%%%%
tmp = spectrum_phi .* U(mod(MU + KK, c * N) + 1);
result = sum(tmp.');

% result=zeros(size(x));
% for kk=1:length(x)
%     for m=0:0
%         PP = mod(mu(kk) + m - K + c*N,c*N)
%         result(kk) = result(kk) + spectrum_phi(kk,m+1)*U(PP+1);
%     end
% end
