%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME DOMAIN INTERPOLATION WITH APPROXIMATE PROLATE SPHEROIDAL WAVEFUNCTION WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = knab2D(ExtrapVal, RHO_IN, COSALPHA_IN, INPUT_DATA, RHO_OUT, COSALPHA_OUT, B1, B2)

% --- Dimensions check
[Nrows, Ncols] = size(INPUT_DATA);
mx = numel(RHO_IN); my = numel(COSALPHA_IN);
if (mx ~= Ncols || my ~= Nrows) && ~isequal(size(RHO_IN),size(COSALPHA_IN),size(INPUT_DATA))
    error('The lengths of the X and Y vectors must match Z.');
end
if (Nrows < 2) || (Ncols < 2)
    error('Z must be at least 2-by-2.');
end
if (Nrows < 2) || (Ncols < 2)
    error('Z must be at least 2-by-2.');
end

% --- Mean sampling steps
T1 = mean(diff(RHO_IN(1, :)));
T2 = mean(diff(COSALPHA_IN(:, 1)));

if T1 > 1/B1
    error('Nyquist frequency is bigger that sampling frequency');    
end
if T2 > 1/B2
    error('Nyquist frequency is bigger that sampling frequency');    
end

P = 18;                                 % --- 2*P+1 is the number of samples involved in the sampling series
Tw1 = P * T1;
Tw2 = P * T2;

result = zeros(size(RHO_OUT, 1), size(COSALPHA_OUT, 2));

for j = 1: size(RHO_OUT, 1)
    for k = 1: size(COSALPHA_OUT, 2)
        
         % --- Time domain approximate prolate spheroidal wave function
         wap1 = sinc((1 / T1 - B1) * sqrt((RHO_OUT(j, k)      - RHO_IN(1, :))    .^2 - Tw1^2)) / sinc(1i * (1 / T1 - B1) * Tw1);
         wap2 = sinc((1 / T2 - B2) * sqrt((COSALPHA_OUT(j, k) - COSALPHA_IN(:,1)).^2 - Tw2^2)) / sinc(1i * (1 / T2 - B2) * Tw2);

         % --- Time domain ectangular window
         g01 = sinc(1 / T1 * (RHO_OUT(j,k)      - RHO_IN(1,:)));
         g02 = sinc(1 / T2 * (COSALPHA_OUT(j,k) - COSALPHA_IN(:,1)));
         
         ke1 = wap1 .* g01;
         ke2 = wap2 .* g02;
         
         [KE1, KE2]     = meshgrid(ke1, ke2);  
         result(j,k)    = result(j,k) + sum(sum(INPUT_DATA .* KE1 .* KE2));
         
    end
end
