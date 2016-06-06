function fval = InvScatt2DNearFieldMultifreqFunc(Coeff_legendre)

global leg leg_der N_sampling_points a NDF_x NDF_z d X_prime Z_prime prolate_spheroidals_x prolate_spheroidals_z beta numFrequencies transmission_reflection

N_unknowns = length(Coeff_legendre);

x       = zeros(1, N_sampling_points / 2);
x_der   = zeros(1, N_sampling_points / 2);
for k = 1 : N_unknowns,
    x       = x     + Coeff_legendre(k) * leg(k, :);
    x_der   = x_der + Coeff_legendre(k) * leg_der(k, :);
end
if (min(x_der) < 0)
    x = x - min(x_der) * leg(1, :);
end
temp = x;
temp = temp - min(temp);
temp = temp(2 : length(temp));
x = [-fliplr(temp) 0 temp];
% x = [-fliplr(x(2 : length(x))) x];
% x = [-fliplr(x(2 : length(x))) 0 x(2 : length(x))];
x = a * x / max(abs(x));

Radiation_Operator = zeros(N_sampling_points * numFrequencies, NDF_x * NDF_z);

for s = 1 : numFrequencies,
    for k = 1 : N_sampling_points - 1,
    
        R = sqrt((x(k) - X_prime).^2 + (d - Z_prime).^2);

        for p = 1 : NDF_x,
            for q = 1 : NDF_z,
                [PROLATE_X, PROLATE_Z] = meshgrid(prolate_spheroidals_x(p, :), prolate_spheroidals_z(q, :).');
                PROLATE_XZ = PROLATE_X .* PROLATE_Z;
                Radiation_Operator((s - 1) * (N_sampling_points - 1) + k, (p - 1) * NDF_z + q) = sum(sum((besselh(0, 2, beta(s) * R) .* exp(transmission_reflection * 1i * beta(s) * Z_prime) .* PROLATE_XZ)));
            end
        end 
    end
end

S = svd(Radiation_Operator);

fval = -sum(S / S(1))
