function result=NFFT1_1D_Optimized(data,x)

N=length(data);

%--- Algorithm parameters
c=2;                                                    % Oversampling factor >=1 (typically 2)
K=3;                                                    % 2K+1 interpolation samples (N should be >> K)
alfa=(2-1/c)*pi-0.01;                                   % Half-size of the support of the interpolation window
alfa_prime=((2-1/c)*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPATIAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk=-N/2:N/2-1;
xi=2*pi*kk/(c*N);

% --- Calculating Prolate Spheroidal Wave Functions (PSWFs), Legendre Polynomials
% and expansion coefficients of Legendre Polynomials yielding the (PSWFs)
% Space_Bandwidth_Product=0.991525423728813*((2*pi-pi/c)*K);
Space_Bandwidth_Product=(1*(2*pi-pi/c)*K);
[PSWFs P V_even V_odd K_Leg] = S0n(Space_Bandwidth_Product, 12, xi / (2 * pi - pi / c), 1e-30);
PSWFs = PSWFs.';

load Result_c_2_NumLegPol_12_K_3_SBPfactor_1.mat

mm=-K:K;

% --- Recovering the "optimal" spatial window by the Legendre Polynomials
phi = zeros(1, length(xi));
for p=0:2:K_Leg-1
    phi = phi + sum(V_even(p/2+1,1:length(c_opt)).*c_opt) * P(p+1,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPECTRAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu=round(c*x);  

[KK,MU]=meshgrid(-K:K,mu);

spectrum_phi = zeros(length(x),length(mm));
[XX,MM]=meshgrid(x,mm);
for p=0:2:K_Leg-1
    temp = sqrt(pi./(2.*alfa_prime*abs(c*XX-(MM+round(c*XX))))) .* besselj(p+0.5,alfa_prime*abs(c*XX-(MM+round(c*XX))));
    indices = abs(c*XX-(MM+round(c*XX))) == 0;
    temp(indices) = sqrt(pi) * (0.5 * alfa_prime*abs(c*XX(indices)-(MM(indices)+round(c*XX(indices))))) .^ p ./ (2 * gamma(p + 1.5));
    spectrum_phi = spectrum_phi + alfa_prime * (sqrt(p+0.5)/(2*pi)) * sum(V_even(p/2+1,1:length(c_opt)).*c_opt) * ((2 * (-1i)^p .* temp)).';
end

%%%%%%%%%%
% STEP 1 %
%%%%%%%%%%

u=zeros(1,c*N);
u((c-1)*N/2+1:(c+1)*N/2)=data./phi;

%%%%%%%%%%
% STEP 2 %
%%%%%%%%%%

U=fft(ifftshift(u));

%%%%%%%%%%
% STEP 3 %
%%%%%%%%%%

tmp=spectrum_phi.*U(mod(MU+KK,c*N)+1);
result=sum(tmp.');

% result=zeros(size(x));
% for kk=1:length(x)
%     for m=0:0
%         PP = mod(mu(kk) + m - K + c*N,c*N)
%         result(kk) = result(kk) + spectrum_phi(kk,m+1)*U(PP+1);
%     end
% end
