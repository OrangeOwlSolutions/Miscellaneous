function [f, xsi_min, xsi_max]= expansiveFuncGen(f0)
   
    % --- Constructing the expansive function f starting from f0
    global C_ord xsi N alfa
    
    Df0 = polyder(f0);
    xsi_refined = linspace(xsi(1), xsi(N), 5 * N);
    func = polyval(Df0, xsi_refined);
    
    [min_Df0, min_ind] = min(func);
    [max_Df0, max_ind] = max(func);
    xsi_min = xsi_refined(min_ind);
    xsi_max = xsi_refined(max_ind);
    
    if abs(max_Df0 - min_Df0) < 1e-3
        f = zeros(1, C_ord);
        f(C_ord - 1) = 1;
    else
        if max_Df0 > alfa
            f = (alfa - 1) / (max_Df0 - min_Df0) * f0;
            f(C_ord - 1) = f(C_ord - 1) - min_Df0 * (alfa - 1) / (max_Df0 - min_Df0) + 1;
        else
            f = f0;
            f(C_ord - 1) = f(C_ord - 1) + (1 - min_Df0);
        end
    end
end
