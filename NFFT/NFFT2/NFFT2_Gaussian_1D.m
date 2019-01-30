function result=NFFT2_Gaussian_1D(data, x, M, c, K)

N = length(data);

%--- Algorithm parameters
alfa=(2-1/c)*pi-0.01;                                   % Half-size of the support of the interpolation window

%--- Algorithm parameters
% b = log(10^(0.135));
% c = 2;                                                  % Oversampling factor >=1 (typically 2)
% K = fix(2 * b * pi);                                    % 2K+1 interpolation samples (N should be >> K) (parameter q/2 in Dutt & Rokhlin)
% b = K / (2 * pi);

b = (2 * c / (2 * c - 1)) * K / (4 * pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPECTRAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[KK, X] = meshgrid(-K : K, x);
M_U     = round(c * X);                                          
phi_cap = exp( -(c * X - (M_U + KK)) .^ 2 / (4 * b)) / (2 * sqrt(b * pi));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPATIAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk = -M / 2 : M / 2 - 1;
xi = 2 * pi * kk / (c * M);
phi = exp(-b * xi .^ 2);

%%%%%%%%%%
% STEP 1 %
%%%%%%%%%%

%--- Old version
% u=zeros(1,c*M);
% PP=mod(M_U+KK+M*c/2,c*M)+1;
% for k=1:N 
%     u(PP(k,(-K:K)+K+1))=u(PP(k,(-K:K)+K+1))+data(k)*phi_cap(k,(-K:K)+K+1);
% end

% --- New version
u = zeros(1, c * M);
for s=-c*M/2:c*M/2-1
    for l=1:N
        phi_cap=0;
        if (abs(round(c*x(l))-s)<=K)
            phi_cap = exp( -(c * x(l) - s) .^ 2 / (4 * b)) / (2 * sqrt(b * pi));
        elseif (abs(round(c*x(l))+c*M-s)<=K)
            phi_cap = exp( -(c * x(l) - s + c * M) .^ 2 / (4 * b)) / (2 * sqrt(b * pi));
        elseif (abs(round(c*x(l))-c*M-s)<=K)
            phi_cap = exp( -(c * x(l) - s - c * M) .^ 2 / (4 * b)) / (2 * sqrt(b * pi));
        end
        u(s+c*M/2+1)=u(s+c*M/2+1)+data(l)*phi_cap;
    end
end

%%%%%%%%%%
% STEP 2 %
%%%%%%%%%%

FFT_Mc = fftshift(fft(ifftshift(u)));

%%%%%%%%%%
% STEP 3 %
%%%%%%%%%%

result = FFT_Mc(1 + (c - 1) * M / 2 : (c + 1) * M / 2) ./ phi;

