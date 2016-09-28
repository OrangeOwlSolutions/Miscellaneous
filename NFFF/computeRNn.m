function RNn = computeRNn(nmax, N, c, r)

% --- For a fixed N, returns RNn for 0 <= n <= nmax - 1

for n = 0 : nmax - 1,
    
    % --- Upper diagonal
    if (n < nmax - 1)
        hNn     = (2 * (2 * (n + 1) + N + 1) * nchoosek((n + 1) + N, (n + 1))^2)^0.5;
        hNn_m1  = (2 * (2 *  n       + N + 1) * nchoosek(n + N,       n     )^2)^0.5;
        gammaNn_m1      = -(hNn / hNn_m1) * (n + 1)^2 /((2 * (n + 1) + N) * (2 * (n + 1) + N + 1));
        B(n + 1, n + 2) = -c^2 * gammaNn_m1;
    end
    
    % --- Main diagonal
    kNn         = (N + 2 * n + 0.5) * (N + 2 * n + 1.5);
    if (N == 0)
        gamma0Nn = 0.5;
    else
        gamma0Nn    = (2 * n * (n + 1) + N * (2 * n + N + 1)) / ((2 * n + N) * (2 * n + N + 2));
    end
    B(n + 1, n + 1) = -(kNn + c^2 * gamma0Nn);
 
    % --- Lower diagonal
    if (n < nmax - 1)
        hNn     = (2 * (2 * n + N + 1) * nchoosek(n + N, n)^2)^0.5;
        hNn_p1  = (2 * (2 * (n + 1) + N + 1) * nchoosek((n + 1) + N, (n + 1))^2)^0.5;
        gammaNn_p1          = -((n + N + 1)^2 / ((2 * n + N + 1) * (2 * n + N + 2))) * hNn / hNn_p1;
        B(n + 2, n + 1) = -c^2 * gammaNn_p1;
    end

end

[D, CHI] = eig(B);

Jacobi_Matrix = Jacobi_PN0n(nmax, N, (1 - 2 * r .^2));

RNn = zeros(nmax, length(r));
for n = 0 : nmax - 1,
    for k = 0 : nmax - 1,
        hNk     = (2 * (2 * k + N + 1) * nchoosek(k + N, k)^2)^0.5;
%         TNk     = hNk * r.^(N + 0.5) .* (nchoosek(k + N, k) ^(-1)) .* Jacobi_Matrix(k + 1, :);
        TNk     = hNk * r.^(N) .* (nchoosek(k + N, k) ^(-1)) .* Jacobi_Matrix(k + 1, :);
        RNn(n + 1, :) = RNn(n + 1, :) + D(k + 1, n + 1) * TNk;
    end
end
    
% --- Ordering 
chi = diag(CHI);
[~, indices] = sort(chi, 'descend');
RNn = RNn(indices, :);
