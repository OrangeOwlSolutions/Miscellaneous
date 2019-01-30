function result=NFFT1_2D(data,x1,x2)

N1=size(data,1);
N2=size(data,2);

M=length(x1);

%% --- Algorithm parameters
c=2;                                                    
K=6;                                                    
alfa=(2-1/c)*pi-0.01;                                   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPATIAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk1=-N1/2:N1/2-1;
kk2=-N2/2:N2/2-1;
xi1=2*pi*kk1/(c*N1);
xi2=2*pi*kk2/(c*N2);
phi1=besseli(0,K*sqrt(alfa.^2-xi1.^2));
phi2=besseli(0,K*sqrt(alfa.^2-xi2.^2));
[PHI2,PHI1]=meshgrid(phi2,phi1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPECTRAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu1=round(c*x1);                                          
mu2=round(c*x2);                                          

[KK1,MU1]=meshgrid(-K:K,mu1);
X1=x1'*ones(1,2*K+1);
P1=sqrt(K.^2-(c*X1-(MU1+KK1)).^2);
spectrum_phi1=(1/pi)*sinh(alfa*P1)./P1;
spectrum_phi1(find(P1==0))=alfa/pi;

[KK2,MU2]=meshgrid(-K:K,mu2);
X2=x2'*ones(1,2*K+1);
P2=sqrt(K.^2-(c*X2-(MU2+KK2)).^2);
spectrum_phi2=(1/pi)*sinh(alfa*P2)./P2;
spectrum_phi2(find(P2==0))=alfa/pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: scaling and zero padding %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u=zeros(c*N1,c*N2);
u((c-1)*N1/2+1:(c+1)*N1/2,(c-1)*N2/2+1:(c+1)*N2/2)=data./(PHI1.*PHI2);

%%%%%%%%%%
% STEP 2 %
%%%%%%%%%%

U=fft2(ifftshift(u));

%%%%%%%%%%
% STEP 3 %
%%%%%%%%%%

result=zeros(1,M);
for l=1:M,
    for m1=-K:K,
        result(l)=result(l)+sum(spectrum_phi1(l,m1+K+1)*spectrum_phi2(l,1:2*K+1).*...
                                                      U(mod(mu1(l)+m1,c*N1)+1,mod(mu2(l)+(-K:K),c*N2)+1));
    end
end

end
