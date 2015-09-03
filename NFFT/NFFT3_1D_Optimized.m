function result = NFFT3_1D_Optimized(data, x, s)

xsi_max = pi;
N = length(data);
M = N;
S = max(abs(s));

%--- Algorithm parameters
c=2;                                                    % Oversampling factor >=1 (typically 2)
K=3;                                                    % 2K+1 interpolation samples (N should be >> K)
alfa=(2-1/c)*pi-0.01;                                   % Half-size of the support of the interpolation window
alfa_prime=((2-1/c)*pi);

x1 = S * x / pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPATIAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% kk=-M/2:M/2-1;
% xi=2*xsi_max*kk/(c*M);

% --- Calculating Prolate Spheroidal Wave Functions (PSWFs), Legendre Polynomials
% and expansion coefficients of Legendre Polynomials yielding the (PSWFs)
% Space_Bandwidth_Product=0.991525423728813*((2*pi-pi/c)*K);
Space_Bandwidth_Product=(1*(2*pi-pi/c)*K);
[PSWFs P V_even V_odd K_Leg] = S0n(Space_Bandwidth_Product, 12, (pi * s / (c * S)) / (2 * pi - pi / c), 1e-30);
PSWFs = PSWFs.';

load Result_c_2_NumLegPol_12_K_3_SBPfactor_1.mat

% --- Recovering the "optimal" spatial window by the Legendre Polynomials
phi2 = zeros(1, length(s));
for p=0:2:K_Leg-1
    phi2 = phi2 + sum(V_even(p/2+1,1:length(c_opt)).*c_opt) * P(p+1,:);
end

%%%%%%%%%%%%%%%
% CONVOLUTION %
%%%%%%%%%%%%%%%

% --- Old version
% u = zeros(1, c * M);
% PP = mod(MU1 + KK + M * c / 2, c * M) + 1;
% for k = 1 : N 
%     u(PP(k, (-K : K) + K + 1)) = u(PP(k, (-K : K) + K + 1)) + data(k) * spectrum_phi1(k, (-K : K) + K + 1);
% end

u=zeros(1,c*M);
for k=-M:M-1
    for l=1:N
        phi_cap(l,k+M+1)=0;
        if (abs(round(c*x1(l))-k)<=K)
            for p=0:2:K_Leg-1
                phi_cap(l,k+M+1) = phi_cap(l,k+M+1) + alfa_prime * (sqrt(p+0.5)/(2*pi)) * sum(V_even(p/2+1,1:length(c_opt)).*c_opt) * ((2 * (-1i)^p .* sqrt(pi./(2.*alfa_prime*abs(c*x(l)*S/pi-k))) .* besselj(p+0.5,alfa_prime*abs(c*x(l)*S/pi-k)))).';
            end
        elseif (abs(round(c*x1(l))+c*M-k)<=K)
            for p=0:2:K_Leg-1
                phi_cap(l,k+M+1) = phi_cap(l,k+M+1) + alfa_prime * (sqrt(p+0.5)/(2*pi)) * sum(V_even(p/2+1,1:length(c_opt)).*c_opt) * ((2 * (-1i)^p .* sqrt(pi./(2.*alfa_prime*abs(c*x(l)*S/pi-k+c*M))) .* besselj(p+0.5,alfa_prime*abs(c*x(l)*S/pi-k+c*M)))).';
            end
        elseif (abs(round(c*x1(l))-c*M-k)<=K)
            for p=0:2:K_Leg-1
                phi_cap(l,k+M+1) = phi_cap(l,k+M+1) + alfa_prime * (sqrt(p+0.5)/(2*pi)) * sum(V_even(p/2+1,1:length(c_opt)).*c_opt) * ((2 * (-1i)^p .* sqrt(pi./(2.*alfa_prime*abs(c*x(l)*S/pi-k-c*M))) .* besselj(p+0.5,alfa_prime*abs(c*x(l)*S/pi-k-c*M)))).';
            end
        end
        u(k+M+1)=u(k+M+1)+data(l)*phi_cap(l,k+M+1);
    end
end

%%%%%%%%%%%%%
% NER NUFFT %
%%%%%%%%%%%%%

result = NFFT1_1D_Optimized(u, s * N / (c * S)) ./ phi2;
