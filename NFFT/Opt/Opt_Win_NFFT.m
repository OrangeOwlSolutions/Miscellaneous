clear all
close all
clc

global mm K EXP_REF PFs xsi xsi_full x indices_xsi XSI MM P K_Leg V_even alfa_prime c

%%%%%%%%%%%%%%
% PARAMETERS %
%%%%%%%%%%%%%%
c = 2;                                                % --- Oversampling factor
alfa_prime = ((2 - 1 / c) * pi);                        % --- Support of Phi

K = 3;                                                  % --- Support of Phi_hat

Max_Num_PFs = 12;                                       % --- Maximum number of PFs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLING IN THE x VARIABLE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 500;                                               % --- Number of samples in the x variable
xM = 10 * pi;                                           % --- Maximum value of the x variable
x = linspace(-xM, xM, Nx);                              % --- x variable samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLING IN THE xsi VARIABLE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nxsi_full = 1800;                                       % --- Number of samples in the xsi variable
xsi_full = linspace(-(2 * pi - 1 * pi / c), 2 * pi - 1 * pi / c, Nxsi_full);
indices_xsi = find(abs(xsi_full) <= pi / c);
xsi = xsi_full(indices_xsi);

mm = -K : K;

%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCE EXPONENTIAL %
%%%%%%%%%%%%%%%%%%%%%%%%%
[XSI2, X]   = meshgrid(xsi, x);
[XSI, MM]   = meshgrid(xsi, mm);
EXP_REF     = exp(-1i * XSI2 .* X);

%%%%%%%
% PFs %
%%%%%%%
% SB_Product=0.991525423728813*((2*pi-pi/c)*K);
SBP_factor                      = 1.;
SB_Product                      = SBP_factor * ((2 * pi - pi / c) * K);
[PFs P V_even V_odd K_Leg]      = S0n(SB_Product, 2 * Max_Num_PFs, xsi_full / (2 * pi - pi / c), 1e-30);
PFs = PFs.';

%%%%%%%%%%%%%%%%
% OPTIMIZATION %
%%%%%%%%%%%%%%%%
options = optimset('Display', 'iter');

% --- Full optimization
c_init = [1 0];
for p=1:Max_Num_PFs,
    
    [c_opt, max_err] = fminsearch(@functional_to_be_optimized, c_init, options); 
    
    c_init = [c_opt 0];
    
end

filename = strcat('Result_c_',num2str(c),'_NumPFs_',num2str(Max_Num_PFs),'_K_',num2str(K),'_SBPfactor_',num2str(SBP_factor),'.mat');
save(filename,'c_opt')

functional_to_be_optimized(c_opt)

%%%%%%%%%%%%%%%%%%%%%%%
% KB WINDOW AND ERROR %
%%%%%%%%%%%%%%%%%%%%%%%
[XXSI, MMM] = meshgrid(xsi, mm);
exp_KB=zeros(length(x),length(xsi));
for k=1:length(x),
    exp_KB(k,:)=exp_KB(k,:)+((2*pi)^(-0.5))*sum(Evaluate_Phi_hat(x(k)-(MMM+round(x(k))),c,K).*...
                                                  exp(-1i*(MMM+round(x(k))).*XXSI))./...
                                                  Evaluate_Phi(xsi,c,K);
end
rms_KB = sum(sum(abs(EXP_REF-exp_KB)))

phi_KB = Evaluate_Phi(xsi,c,K);
phi_KB2 = Evaluate_Phi(xsi_full,c,K);

%%%%%%%%%%%%%%%%%%%%
% OPTIMIZED WINDOW %
%%%%%%%%%%%%%%%%%%%%
% --- Recovering the "optimal" spatial window 
Phi = zeros(length(xsi_full),1);
for p=0:2:K_Leg-1
    Phi = Phi + sum(V_even(p/2+1,1:length(c_opt)).*c_opt) * P(p+1,:).';
end
Phi2=Phi(indices_xsi);

%%%%%%%%%%
% GRAPHS %
%%%%%%%%%%
% —- For drawings only
set(gca,'FontSize',14,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')

figure(1)
plot(xsi / (pi / c), abs(phi_KB) / max(abs(phi_KB)), 'r--', 'LineWidth', 2)
hold on
plot(xsi / (pi / c), abs(Phi2) / max(abs(Phi2)), 'LineWidth', 2)
hold off
xlabel('\xi / (\pi/c)')
ylabel('\phi(\xi)')

figure(2)
plot(xsi_full / (2 * pi - 1 * pi / c), abs(phi_KB2) / max(abs(phi_KB2)), 'r--', 'LineWidth', 2)
hold on
plot(xsi_full / (2 * pi - 1 * pi / c), abs(Phi) / max(abs(Phi)), 'LineWidth', 2)
hold off
xlabel('\xi / (2\pi-\pi/c)')
ylabel('\phi(\xi)')


