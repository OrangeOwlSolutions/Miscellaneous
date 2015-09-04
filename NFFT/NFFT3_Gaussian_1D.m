function result = NFFT3_Gaussian_1D(data, x, s)

N = length(data);
M = N;
S = max(abs(s));

%--- Algorithm parameters
b = log(10^(0.135));
c = 2;                                                  % Oversampling factor >=1 (typically 2)
K = fix(2 * b * pi);                                    % 2K+1 interpolation samples (N should be >> K) (parameter q/2 in Dutt & Rokhlin)

x1 = S * x / pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPATIAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi2 = exp(-b * (pi * s / (c * S)) .^2);

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
        phi_cap=0;
        if (abs(round(c*x1(l))-k)<=K)
            phi_cap = exp( -(c * x1(l) - k) .^ 2 / (4 * b)) / (2 * sqrt(b * pi));
        elseif (abs(round(c*x1(l))+c*M-k)<=K)
            phi_cap = exp( -(c * x1(l) - k + c * M) .^ 2 / (4 * b)) / (2 * sqrt(b * pi));
        elseif (abs(round(c*x1(l))-c*M-k)<=K)
            phi_cap = exp( -(c * x1(l) - k - c * M) .^ 2 / (4 * b)) / (2 * sqrt(b * pi));
        end
        u(k+M+1)=u(k+M+1)+data(l)*phi_cap;
    end
end

%%%%%%%%%%%%%
% NER NUFFT %
%%%%%%%%%%%%%

result = NFFT1_Gaussian_1D(u, s * N / (c * S)) ./ phi2;
