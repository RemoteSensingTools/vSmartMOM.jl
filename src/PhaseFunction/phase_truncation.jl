# Just prototype code right now, will add real functions and AbstractTruncationType

n = 20
μ_max = 730
y = bulk_f₁₁[1:μ_max]
K = transpose(P[1:n,1:μ_max])
Se = Diagonal((ones(length(y)).*y).^2);
β_trun = inv(K'inv(Se)*K)K'inv(Se)*y