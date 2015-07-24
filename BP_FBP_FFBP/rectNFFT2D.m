function result = rectNFFT2D(ExtrapVal, RHO_IN, COSALPHA_IN, INPUT_DATA, RHO_OUT, COSALPHA_OUT, B1, B2)

% --- ExtrapVal:                            Value to be assigned to all output points for which extrapolation is needed                         

% --- Mean sampling steps
T1 = mean(diff(RHO_IN(1, :)));
T2 = mean(diff(COSALPHA_IN(:, 1)));

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

% --- Shifting the output sampling locations
RHO_OUT_SHIFTED         = RHO_OUT - RHO_IN(1);
COSALPHA_OUT_SHIFTED    = COSALPHA_OUT - COSALPHA_IN(1);

len1 = size(RHO_OUT_SHIFTED, 1);
len2 = size(RHO_OUT_SHIFTED, 2);

min_RHO_OUT_SHIFTED         = (min(min(RHO_OUT_SHIFTED)));
min_COSALPHA_OUT_SHIFTED    = (min(min(COSALPHA_OUT_SHIFTED)));
max_RHO_OUT_SHIFTED         = (max(max(RHO_OUT_SHIFTED)));
max_COSALPHA_OUT_SHIFTED    = (max(max(COSALPHA_OUT_SHIFTED)));

N = ((max_RHO_OUT_SHIFTED - min_RHO_OUT_SHIFTED) / (T1 * (1 - 1 / len2)));
M = ((max_COSALPHA_OUT_SHIFTED - min_COSALPHA_OUT_SHIFTED) / (T2 * (1 - 1 / len1)));

if (N<1) || (isnan(N))
    N=1;
end
if (M<1) || (isnan(M))
    M=1;
end

% --- Oversampling detection
a = round(len2 / N);
b = round(len1 / M);

if a < 2
    a = 2;
    ii = round(a * N) - size(RHO_OUT_SHIFTED, 2);
    RHO_OUT_SHIFTED         = [RHO_OUT_SHIFTED,         zeros(size(RHO_OUT_SHIFTED, 1),         ii)];
    COSALPHA_OUT_SHIFTED    = [COSALPHA_OUT_SHIFTED,    zeros(size(COSALPHA_OUT_SHIFTED, 1),    ii)];
end

if b < 2
    b = 2;
    ii = round(b * M) - size(RHO_OUT_SHIFTED, 1);
    RHO_OUT_SHIFTED         = [RHO_OUT_SHIFTED;         zeros(ii, size(RHO_OUT_SHIFTED, 2))];
    COSALPHA_OUT_SHIFTED    = [COSALPHA_OUT_SHIFTED;    zeros(ii, size(COSALPHA_OUT_SHIFTED, 2))];
end

dfn = 1 / (N * T1);
kgn = floor((1 / T1 - B1 / 2) / dfn); 
dfm = 1 / (M * T2);
kgm = floor((1 / T2 - B2 / 2) / dfm); 

%--- NUFFT based Selva's approach to bandlimited signal interpolation with
%rectangular window
result = reshape(NFFT_Based_Selvas_Approach_for_Signal_Interp_2D(INPUT_DATA.' / size(INPUT_DATA.', 1) / size(INPUT_DATA.', 2), T1, T2, a, b, 1, kgn, kgm, RHO_OUT_SHIFTED, COSALPHA_OUT_SHIFTED), size(RHO_OUT_SHIFTED, 1), size(COSALPHA_OUT_SHIFTED, 2));
result = result(1:len1, 1:len2);

% --- Setting ExtrapVal to all extrapolation points
xout = (RHO_OUT(1,:)      < min(RHO_IN(1,:)))      | (RHO_OUT(1,:)      > max(RHO_IN(1,:)));
yout = (COSALPHA_OUT(:,1) < min(COSALPHA_IN(:,1))) | (COSALPHA_OUT(:,1) > max(COSALPHA_IN(:,1)));
result(:,xout) = ExtrapVal;
result(yout,:) = ExtrapVal;

