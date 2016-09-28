function Spectrum = Ea2S(uu, vv, Ea, x, y, delta_x_1, delta_y_1)   

[YY, VV] = meshgrid(y, vv);

[DELTAy] = meshgrid(delta_y_1, vv);            

W1       = exp(1i * (YY .* VV)) .* DELTAy;          

[XX, UU] = meshgrid(x, uu);               

[DELTAx] = meshgrid(delta_x_1, uu);

W2       = exp(1i * (XX .* UU)) .* DELTAx;

Spectrum = W1 * (W2 * Ea')';               
