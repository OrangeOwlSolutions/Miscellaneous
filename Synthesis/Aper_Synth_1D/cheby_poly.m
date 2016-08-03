function chebychev_poly = cheb_poly(C_ord)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Calculates the Chebyshev polynomials up to order C_ord        %
%                      WARNING!                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    chebychev_poly = zeros(C_ord, C_ord);                   % --- Initialize the polynomial coefficient matrix

    chebychev_poly(1, C_ord)     = 1;                       % --- 0-th order polynomial (costant)
    chebychev_poly(2, C_ord - 1) = 1;                       % --- 1-st order polynomial (linear)

    for kk = 3 : C_ord   
        chebychev_poly(kk, :) = 2 * circshift(chebychev_poly(kk - 1, :).', -1).' - chebychev_poly(kk - 2, :);
    end

end
   
