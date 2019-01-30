function [PFs P V_even V_odd K_Leg] = S0n(c, N, x, epsilon)

% --- Chooses the maximum number of polynomials used to represent the PFs

if N<ceil(2*c/pi)
    K_Leg=ceil(2*(ceil(exp(1)*c)+1)+log2(1/epsilon));
else
    lambda_m=exp((N+0.5)*(log(exp(1)*c/4)-log(N+0.5)));
    K_Leg=ceil(2*(ceil(exp(1)*c)+1)+log2(1/epsilon)+log2(1/lambda_m));
end

% --- Generates Lag polynomials and normalizes them

P=Legendre_Polynomials(K_Leg+1,x);

for k=0:K_Leg+1,
    P(k+1,:)=P(k+1,:)*sqrt(k+0.5);
end

% --- Generates the matrix A and its even and odd parts. Evaluates their eigenvectors

A=zeros(K_Leg+2,K_Leg+2);
for k=0:K_Leg-1,
    A(k+1,k+1)=k*(k+1)+((2*k*(k+1)-1)/((2*k+3)*(2*k-1)))*c^2;
    A(k+1,k+2+1)=(((k+2)*(k+1)/((2*k+3)*sqrt((2*k+1)*(2*k+5)))))*c^2;
    A(k+2+1,k+1)=(((k+2)*(k+1))/((2*k+3)*sqrt((2*k+1)*(2*k+5))))*c^2;
end

A=A(1:K_Leg,1:K_Leg);

A_even=A(1:2:size(A,1),1:2:size(A,2));
A_odd=A(2:2:size(A,1),2:2:size(A,2));

[V_even,D_even]=eig(A_even);
[V_odd,D_odd]=eig(A_odd);

% --- Forms the PFs

PFs=zeros(N,length(x));

for m=1:2:N,
    for n=1:size(V_even,1),
        PFs(m,:)=PFs(m,:)+V_even(n,(m-1)/2+1)*P(2*n-1,:);
    end
end

for m=2:2:N,
    for n=1:size(V_odd,1),
        PFs(m,:)=PFs(m,:)+V_odd(n,m/2)*P(2*n,:);
    end
end

