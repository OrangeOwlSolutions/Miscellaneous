function S = batchedMatlabSVD(A, Nrows, Ncols, numMatrices)
  
B = zeros(Nrows, Ncols);
   
for p = 0 : numMatrices - 1
       
    for k = 1 : Nrows * Ncols
        B(k) = A(k + p * Nrows * Ncols);
    end
    
    s = svd(B);
       
    for h = 1 : Ncols
        S(h + p * Ncols) = s(h);
        S = S';
    end
    
end
