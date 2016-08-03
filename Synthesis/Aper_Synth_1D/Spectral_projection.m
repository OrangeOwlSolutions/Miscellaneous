function Projected_spectrum=Spectral_projection(Spectrum)

global Upper_mask Lower_mask

Projected_spectrum=Spectrum;

indices_upper=find(abs(Spectrum) > Upper_mask);
indices_lower=find(abs(Spectrum) < Lower_mask);

Projected_spectrum(indices_upper) = Upper_mask(indices_upper) .* exp(1i * angle(Spectrum(indices_upper)));
Projected_spectrum(indices_lower) = Lower_mask(indices_lower) .* exp(1i * angle(Spectrum(indices_lower)));

