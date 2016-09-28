%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                                                             %
% Calculates Jacobi polynomials P_n^(N,0)(z) for 0 <= n <= nmax               %                    
%                                                                             %
% Uses Rodrigues' formula                                                     % 
%                                                                             %
% P_n^(N,0)(z) = ((-1)^n/(2^n n!)) (1-z)^(-N) d^n/dz^n [(1-z)^N (1-z^2)^n]    %  
%                                                                             %
% checked against http://keisan.casio.com/has10/SpecExec.cgi                  %                
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function Jacobi_Matrix=Jacobi_PN0n(nmax,N,x);

Jacobi_Matrix=zeros(nmax+1,length(x));

Jacobi_Matrix(1,:)=ones(size(x));

if (N>0 && nmax>0)
    for n=1:nmax
        Factor=2^(-n)/factorial(n);
        RootsP1=ones(1,n);
        RootsP2=-ones(1,n);
        RootsP3=ones(1,N);
        CoeffP1=poly(RootsP1);
        CoeffP2=poly(RootsP2);
        CoeffP3=poly(RootsP3);
        CoeffP1P2=conv(CoeffP1,CoeffP2);
        CoeffDer=polyder(CoeffP1P2,CoeffP3);
        for k=1:n-1,
            CoeffDer=polyder(CoeffDer);
        end
        RootsP4=ones(1,N);
        CoeffP4=poly(RootsP4);
        CoeffP5=deconv(CoeffDer,CoeffP4);
        Jacobi_Matrix(n+1,:)=Factor*polyval(CoeffP5,x);
    end
end

if (N==0 && nmax>0)
    for n=1:nmax
        Factor=2^(-n)/factorial(n);
        RootsP1=ones(1,n);
        RootsP2=-ones(1,n);
        CoeffP1=poly(RootsP1);
        CoeffP2=poly(RootsP2);
        CoeffDer=polyder(CoeffP1,CoeffP2);                             % --- First order derivative
        for k=1:n-1,                                                    
           CoeffDer=polyder(CoeffDer);                                 % --- Higher order derivatives
        end
        Jacobi_Matrix(n+1,:)=Factor*polyval(CoeffDer,x);
    end
end

