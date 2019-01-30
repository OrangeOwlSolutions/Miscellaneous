function result = NFFT1_1D(data, x, c, K)

N=length(data);

%--- Algorithm parameters
alfa=(2-1/c)*pi-0.01;                                   % Half-size of the support of the interpolation window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPATIAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk=-N/2:N/2-1;
xi=2*pi*kk/(c*N);
phi=besseli(0,K*sqrt(alfa.^2-xi.^2))/besseli(0,K*alfa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPECTRAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu=round(c*x);  

[KK,MU]=meshgrid(-K:K,mu);
X=x'*ones(1,2*K+1);
P=sqrt(K.^2-(c*X-(MU+KK)).^2);

spectrum_phi=(1/pi)*sinh(alfa*P)./P;
spectrum_phi((P==0))=alfa/pi;
spectrum_phi=spectrum_phi/besseli(0,K*alfa);

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
