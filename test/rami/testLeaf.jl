# Quadrature points in μ
n_μ  = 50;
n_ϕ  = 30;
n_ϕᴸ = 55;
n_μᴸ = 40;

R = FT(0.5)
T = FT(0.5)

FT = Float32
μ,w         = CanopyOptics.gauleg(n_μ, FT(0.0),FT(1.0));
dϕ,  w_azi  = CanopyOptics.gauleg(n_ϕ, FT(0),FT(π));
dϕᴸ, w_aziᴸ = CanopyOptics.gauleg(n_ϕᴸ,FT(0),FT(2π));
θᴸ,wᴸ       = CanopyOptics.gauleg(n_μᴸ,FT(0.0),FT(π/2));
μᴸ = cos.(θᴸ)

# Reshape stuff:

w2 = reshape(w2,1,1,1,80,1);


μⁱⁿ   = reshape(CuArray(μ),  n_μ, 1,   1,  1,   1   );
μᵒᵘᵗ  = reshape(CuArray(-μ),  1,   n_μ, 1,  1,   1   );
_dϕ   = reshape(CuArray(dϕ), 1   ,1,   n_ϕ,1,   1   );
_μᴸ   = reshape(CuArray(μᴸ), 1   ,1,   1,  n_μᴸ,1   );
_dϕᴸ  = reshape(CuArray(dϕᴸ),1  , 1,   1,  1,   n_ϕᴸ);

# Quadrature points:
_w_azi  = reshape(CuArray(w_azi),  1   ,1,   n_ϕ,1,     1   );
wᴸ  = wᴸ .*  sin.(θᴸ)
_wᴸ     = reshape(CuArray(wᴸ),     1   ,1,     1,  n_μᴸ,  1   );
_w_aziᴸ = reshape(CuArray(w_aziᴸ), 1   ,1, 1,  1,       n_ϕᴸ);

 







b = sum(sum(sum(integrand.*_w_aziᴸ, dims=5).*_wᴸ,dims=4).*_w_azi, dims=3);
b2 = sum(sum(integrand.*_w_aziᴸ, dims=5).*_wᴸ,dims=4);
b3 = sum(integrand.*_w_aziᴸ, dims=5);

temp = zeros(size(integrand)[1:3])
for i in CartesianIndices(temp)
    temp[i[1], i[2], i[3]] = (wᴸ' * (Array(integrand[i[1], i[2], i[3],:,:]) * w_aziᴸ));
end
temp = (wᴸ' * (Array(integrand[1,1,1,:,:]) * w_aziᴸ))


b = reshape(b,n_μ,n_μ)


μ₁ = reshape(μ₁, 80,1,1,1,1)


