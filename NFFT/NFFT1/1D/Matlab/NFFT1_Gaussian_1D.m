function result=NFFT1_Gaussian_1D(data, x, c, K)

N=length(data);

%--- Algorithm parameters
% b = log(10^(0.135));
% K = fix(2 * b * pi);                                    % 2K+1 interpolation samples (N should be >> K) (parameter q/2 in D&R page 1380)
% b = K / (2 * pi);

b = (2 * c / (2 * c - 1)) * K / (4 * pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPATIAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk = -N / 2 : N / 2 - 1;
xi = 2 * pi * kk / (c * N);
phi = exp(-b * xi .^ 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPECTRAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = round(c * x);  

[KK, MU] = meshgrid(-K : K, mu);
X = x' * ones(1, 2 * K + 1);
spectrum_phi = exp( -(c * X - (MU + KK)) .^ 2 / (4 * b)) / (2 * sqrt(b * pi));

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
