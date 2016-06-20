function max_err = functional_to_be_optimized(c_in)

global EXP_REF V_even x xsi xsi_full mm indices_xsi MM P K_Leg alfa_prime

% --- K_Leg         : Maximum number of polynomials used to represent the PFs

% --- Recovering the "optimal" spatial window 
Phi = zeros(length(xsi_full), 1);
for p = 0 : 2 : K_Leg-1
    Phi = Phi + sum(V_even(p / 2 + 1, 1 : length(c_in)) .* c_in) * P(p + 1, :).';
end 
Phi2 = Phi(indices_xsi);

PHI = Phi2 * ones(size(x));

PHI_HAT     = zeros(length(x), length(mm));
[XX, MM]    = meshgrid(x, mm);
exp_approx  = zeros(length(xsi), length(x));
for p = 0 : 2 : K_Leg - 1
    % --- Using the Fourier transform of Leg polynomials
    PHI_HAT = PHI_HAT + alfa_prime * (sqrt(p + 0.5) / (2 * pi)) * sum(...
                        V_even(p / 2 + 1, 1 : length(c_in)) .* c_in) * ((2 * (-1i)^p .* sqrt(pi ./ (2 .* alfa_prime * abs(XX - (MM + round(XX))))) .* ...
                        besselj(p + 0.5, alfa_prime * abs(XX - (MM + round(XX)))))).';
end

for k = 1 : length(xsi),
    exp_approx(k,:) = sum((PHI_HAT.') .* exp(-1i * (MM + round(XX)) .* xsi(k)));
end
max_err = sum(sum(abs(EXP_REF - exp_approx.' ./ PHI.')));


