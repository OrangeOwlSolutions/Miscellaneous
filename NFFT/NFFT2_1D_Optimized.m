function result = NFFT2_1D_Optimized(data, x, M)
   
xsi_max=pi;
N=length(data);                                        % Number of input elements
    
%--- Algorithm parameters
c=2;                                                    % Oversampling factor >=1 (typically 2)
K=3;                                                    % 2K+1 interpolation samples (N should be >> K)
alfa=(2-1/c)*pi-0.01;                                   % Half-size of the support of the interpolation window
alfa_prime=((2-1/c)*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPATIAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk=-M/2:M/2-1;
xi=2*xsi_max*kk/(c*M);

% --- Calculating Prolate Spheroidal Wave Functions (PSWFs), Legendre Polynomials
% and expansion coefficients of Legendre Polynomials yielding the (PSWFs)
% Space_Bandwidth_Product=0.991525423728813*((2*pi-pi/c)*K);
Space_Bandwidth_Product=(1*(2*pi-pi/c)*K);
[PSWFs P V_even V_odd K_Leg] = S0n(Space_Bandwidth_Product, 12, xi / (2 * pi - pi / c), 1e-30);
PSWFs = PSWFs.';

load Result_c_2_NumLegPol_12_K_3_SBPfactor_1.mat

% --- Recovering the "optimal" spatial window by the Legendre Polynomials
phi = zeros(1, length(xi));
for p=0:2:K_Leg-1
    phi = phi + sum(V_even(p/2+1,1:length(c_opt)).*c_opt) * P(p+1,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPECTRAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mm=-K:K;

phi_cap = zeros(length(x),length(mm));
[XX,MM]=meshgrid(x,mm);
for p=0:2:K_Leg-1
    temp = sqrt(pi./(2.*alfa_prime*abs(c*XX-(MM+round(c*XX))))) .* besselj(p+0.5,alfa_prime*abs(c*XX-(MM+round(c*XX))));
    indices = abs(c*XX-(MM+round(c*XX))) == 0;
    temp(indices) = sqrt(pi) * (0.5 * alfa_prime*abs(c*XX(indices)-(MM(indices)+round(c*XX(indices))))) .^ p ./ (2 * gamma(p + 1.5));
    phi_cap = phi_cap + alfa_prime * (sqrt(p+0.5)/(2*pi)) * sum(V_even(p/2+1,1:length(c_opt)).*c_opt) * ((2 * (-1i)^p .* temp)).';
end

%%%%%%%%%%
% STEP 1 %
%%%%%%%%%%

% --- Old version
% u=zeros(1,c*M);
% [KK,X]=meshgrid(-K:K,x);
% M_U=round(c*X*xsi_max/pi);                                          
% PP=mod(M_U+KK+M*c/2,c*M)+1;
% for k=1:N 
%     u(PP(k,(-K:K)+K+1))=u(PP(k,(-K:K)+K+1))+data(k)*phi_cap(k,(-K:K)+K+1);
% end

% --- New version
u=zeros(1,c*M);
for s=-M:M-1
    for l=1:N
        phi_cap=0;
        if (abs(round(c*x(l))-s)<=K)
            for p=0:2:K_Leg-1
                phi_cap = phi_cap + alfa_prime * (sqrt(p+0.5)/(2*pi)) * sum(V_even(p/2+1,1:length(c_opt)).*c_opt) * ((2 * (-1i)^p .* sqrt(pi./(2.*alfa_prime*abs(c*x(l)-s))) .* besselj(p+0.5,alfa_prime*abs(c*x(l)-s)))).';
            end
        elseif (abs(round(c*x(l))+c*M-s)<=K)
            for p=0:2:K_Leg-1
                phi_cap = phi_cap + alfa_prime * (sqrt(p+0.5)/(2*pi)) * sum(V_even(p/2+1,1:length(c_opt)).*c_opt) * ((2 * (-1i)^p .* sqrt(pi./(2.*alfa_prime*abs(c*x(l)-s+c*M))) .* besselj(p+0.5,alfa_prime*abs(c*x(l)-s+c*M)))).';
            end
        elseif (abs(round(c*x(l))-c*M-s)<=K)
            for p=0:2:K_Leg-1
                phi_cap = phi_cap + alfa_prime * (sqrt(p+0.5)/(2*pi)) * sum(V_even(p/2+1,1:length(c_opt)).*c_opt) * ((2 * (-1i)^p .* sqrt(pi./(2.*alfa_prime*abs(c*x(l)-s-c*M))) .* besselj(p+0.5,alfa_prime*abs(c*x(l)-s-c*M)))).';
            end
        end
        u(s+M+1)=u(s+M+1)+data(l)*phi_cap;
    end
end

%%%%%%%%%%
% STEP 2 %
%%%%%%%%%%

% FFT_Mc=pi/xsi_max*M*c*ifftshift(ifft(fftshift(u)));
FFT_Mc=pi/xsi_max*fftshift(fft(ifftshift(u)));

%%%%%%%%%%
% STEP 3 %
%%%%%%%%%%

result=FFT_Mc(1+(c-1)*M/2:(c+1)*M/2)./phi;

end 
