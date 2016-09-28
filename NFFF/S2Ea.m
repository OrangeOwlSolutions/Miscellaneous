function Ea = S2Ea(uu, vv, Spectrum, x, y, delta_uu, delta_vv)                    

[VV, YY] = meshgrid(vv, y);
                                  
[DELTAv] = meshgrid(delta_vv, y);      

W1       = exp(-1i * (YY .* VV)) .* DELTAv;      

[UU, XX] = meshgrid(uu, x);                  
                                  
[DELTAu] = meshgrid(delta_uu, x);      

W2       = exp(-1i * (XX .* UU)) .* DELTAu;       

Ea       = (1 / (4 * pi^2)) * (W1 * (W2 * Spectrum')'); 
