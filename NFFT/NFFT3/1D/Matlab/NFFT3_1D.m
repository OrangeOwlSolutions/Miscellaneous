function result = NFFT3_1D(data, x, s, c, K)

N = length(data);
M = N;
S = max(abs(s));

%--- Algorithm parameters
alfa=(2-1/c)*pi-0.01;                                   % Half-size of the support of the interpolation window

x1 = S * x / pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPATIAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi2 = besseli(0, K * sqrt(alfa .^2 - (pi * s / (c * S)) .^2 ));

%%%%%%%%%%%%%%%
% CONVOLUTION %
%%%%%%%%%%%%%%%

% --- Old version
% u = zeros(1, c * M);
% PP = mod(MU1 + KK + M * c / 2, c * M) + 1;
% for k = 1 : N 
%     u(PP(k, (-K : K) + K + 1)) = u(PP(k, (-K : K) + K + 1)) + data(k) * spectrum_phi1(k, (-K : K) + K + 1);
% end

u = zeros(1, 2 * M);
phi_cap = zeros(N, 2 * M);
for k = -M : M - 1
    for l = 1 : N
        phi_cap(l, k+M+1) =0;
        if (abs(round(c*x1(l))-k)<=K)
            pp=sqrt(K.^2-(c*x1(l)-k).^2);
            phi_cap(l, k+M+1)=(1/pi)*sinh(alfa*pp)./pp;
            if pp==0 phi_cap(l, k+M+1)=alfa/pi; end
%         elseif (abs(round(c*x1(l))+2*M-k)<=K)
%             pp=sqrt(K.^2-(c*x1(l)-k+2*M).^2);
%             phi_cap=(1/pi)*sinh(alfa*pp)./pp;
%             if pp==0 phi_cap=alfa/pi; end
%         elseif (abs(round(c*x1(l))-2*M-k)<=K)
%             pp=sqrt(K.^2-(c*x1(l)-k-2*M).^2);
%             phi_cap=(1/pi)*sinh(alfa*pp)./pp;
%             if pp==0 phi_cap=alfa/pi; end
        end
%         u(k+M/2+1)=u(k+M/2+1)+data(l)*phi_cap / besseli(0, K * alfa);
        u(k+M+1)=u(k+M+1)+data(l)*phi_cap(l, k+M+1);
    end
end

% u = zeros(1, c * M);
% phi_cap = zeros(N, c * M);
% for k = -c*M/2 : c*M/2 - 1
%     for l = 1 : N
%         phi_cap(l, k+c*M/2+1) =0;
%         if (abs(round(c*x1(l))-k)<=K)
%             pp=sqrt(K.^2-(c*x1(l)-k).^2);
%             phi_cap(l, k+M+1)=(1/pi)*sinh(alfa*pp)./pp;
%             if pp==0 phi_cap(l, k+c*M/2+1)=alfa/pi; end
%         elseif (abs(round(c*x1(l))+c*M-k)<=K)
%             pp=sqrt(K.^2-(c*x1(l)-k+c*M).^2);
%             phi_cap=(1/pi)*sinh(alfa*pp)./pp;
%             if pp==0 phi_cap=alfa/pi; end
%         elseif (abs(round(c*x1(l))-c*M-k)<=K)
%             pp=sqrt(K.^2-(c*x1(l)-k-c*M).^2);
%             phi_cap=(1/pi)*sinh(alfa*pp)./pp;
%             if pp==0 phi_cap=alfa/pi; end
%         end
%         u(k+c*M/2+1)=u(k+c*M/2+1)+data(l)*phi_cap(l, k+c*M/2+1);
%     end
% end

%%%%%%%%%%%%%
% NER NUFFT %
%%%%%%%%%%%%%

result = NFFT1_1D(u, s * N / (c * S), c, K) ./ phi2;

