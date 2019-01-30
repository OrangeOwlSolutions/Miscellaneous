function result = NFFT2_1D_Opt(data, x, M, c, K, Max_Num_PFs, SBP_factor)
   
xsi_max = pi;
N = length(data);                                       % Number of input elements
    
%--- Algorithm parameters
alfa_prime = ((2 - 1 / c) * pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPATIAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk = -M / 2 : M / 2 - 1;
xi = 2 * xsi_max * kk / (c * M);

%%%%%%%
% PFs %
%%%%%%%
SBP_product=(SBP_factor * (2 * pi - pi / c) * K);
[PSWFs P V_even V_odd K_Leg] = S0n(SBP_product, 2 * Max_Num_PFs, xi / (2 * pi - pi / c), 1e-30);
PSWFs = PSWFs.';

filename = strcat('Result_c_', num2str(c), '_NumPFs_', num2str(Max_Num_PFs), '_K_',num2str(K), '_SBPfactor_', num2str(SBP_factor), '.mat');
load(filename)

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
u=zeros(1, c * M);
for s=-c*M/2:c*M/2-1
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
        u(s+c*M/2+1)=u(s+c*M/2+1)+data(l)*phi_cap;
    end
end

% u = zeros(1, c * M);
% for k=-c*M/2:c*M/2-1
%     for l=1:N
%         phi_cap=0;
%         if (abs(round(c*x(l))-k)<=K)
%             pp=sqrt(K.^2-(c*x(l)-k).^2);
%             phi_cap=(1/pi)*sinh(alfa*pp)./pp;
%             if pp==0 phi_cap=alfa/pi; end
%         elseif (abs(round(c*x(l))+c*M-k)<=K)
%             pp=sqrt(K.^2-(c*x(l)-k+c*M).^2);
%             phi_cap=(1/pi)*sinh(alfa*pp)./pp;
%             if pp==0 phi_cap=alfa/pi; end
%         elseif (abs(round(c*x(l))-c*M-k)<=K)
%             pp=sqrt(K.^2-(c*x(l)-k-c*M).^2);
%             phi_cap=(1/pi)*sinh(alfa*pp)./pp;
%             if pp==0 phi_cap=alfa/pi; end
%         end
%         u(k+c*M/2+1)=u(k+c*M/2+1)+data(l)*phi_cap;
%     end
% end

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
