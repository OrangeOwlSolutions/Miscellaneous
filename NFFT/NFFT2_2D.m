function transf = NFFT2_2D(data, x1, x2, M1, M2)
   
N=length(data);                                 % Number of input elements                              

% --- Algorithm parameters
c=2;                                            % Oversampling factor >=1 (typically 2)    
K=6;                                            % 2K+1 interpolation samples (N should be >> K)
alfa=(2-1/c)*pi-0.01;                           % Half-size of the support of the interpolation window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPATIAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu1=round(c*x1);                                          
mu2=round(c*x2);                                         
   
[KK1,MU1]=meshgrid(-K:K,mu1);
X1=x1'*ones(1,2*K+1);
P1=sqrt(K.^2-(c*X1-(MU1+KK1)).^2);
phi_cap1=(1/pi)*sinh(alfa*P1)./P1;
phi_cap1(find(P1==0))=alfa/pi;

[KK2,MU2]=meshgrid(-K:K,mu2);
X2=x2'*ones(1,2*K+1);
P2=sqrt(K.^2-(c*X2-(MU2+KK2)).^2);
phi_cap2=(1/pi)*sinh(alfa*P2)./P2;
phi_cap2(find(P2==0))=alfa/pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE SPECTRAL WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk1=-M1/2:M1/2-1;
xi1=2*pi*kk1/(c*M1);
phi1=besseli(0,K*sqrt(alfa.^2-xi1.^2));
    
kk2=-M2/2:M2/2-1;
xi2=2*pi*kk2/(c*M2);
phi2=besseli(0,K*sqrt(alfa.^2-xi2.^2));

[PHI2,PHI1]=meshgrid(phi2,phi1);

%%%%%%%%%%
% STEP 1 %
%%%%%%%%%%

u=zeros(c*M1,c*M2);
for k=1:N
    for l1=-K:K
        for l2=-K:K,
            p1=mod(mu1(k)+l1+M1*c/2,c*M1)+1;
            p2=mod(mu2(k)+l2+M2*c/2,c*M2)+1;
            u(p1,p2)=u(p1,p2)+data(k)*phi_cap1(k,l1+K+1)*phi_cap2(k,l2+K+1);
        end
    end
end

%%%%%%%%%%
% STEP 2 %
%%%%%%%%%%

% FFT_Mc=M1*M2*(c*c)*ifftshift(ifft2(fftshift(u)));
FFT_Mc=fftshift(fft2(ifftshift(u)));
% FFT_Mc=fft2(u);

%%%%%%%%%%
% STEP 3 %
%%%%%%%%%%

transf=FFT_Mc(1+(c-1)*M1/2:(c+1)*M1/2,1+(c-1)*M2/2:(c+1)*M2/2)./(PHI1.*PHI2);

end 
