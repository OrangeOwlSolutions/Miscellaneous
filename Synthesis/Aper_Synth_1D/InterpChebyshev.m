function InterChebyPoly = InterpChebyshev(c)

    global N_nodes C_ord
    
    ind_j =  0 : (C_ord - 1);  
    ind_k = (1 : N_nodes) - 1 / 2;
    [mat_j, mat_k] = meshgrid(ind_j, ind_k);
    
    COS = cos(pi * (mat_j .* mat_k / N_nodes));
    coeff = (2 / N_nodes) * c * COS;        
    chebyPoly = cheby_poly(C_ord);
    
    InterChebyPoly = coeff * chebyPoly;
    InterChebyPoly(C_ord) = InterChebyPoly(C_ord) - .5 * coeff(1);          
   
end
