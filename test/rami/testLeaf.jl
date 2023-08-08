using CanopyOptics
using CUDA
FT = Float64
# Quadrature points in μ
n_μ  = 50;
n_ϕ  = 30;
n_ϕᴸ = 55;
n_μᴸ = 40;

R = FT(0.5)
T = FT(0.5)


μ,w         = CanopyOptics.gauleg(n_μ, FT(0.0),FT(1.0));
dϕ,  w_azi  = CanopyOptics.gauleg(n_ϕ, FT(0),FT(π));
dϕᴸ, w_aziᴸ = CanopyOptics.gauleg(n_ϕᴸ,FT(0),FT(2π));
θᴸ,wᴸ       = CanopyOptics.gauleg(n_μᴸ,FT(0.0),FT(π/2));
μᴸ = cos.(θᴸ)

# Reshape stuff:
μⁱⁿ   = reshape(CuArray(μ),  n_μ, 1,   1,  1,   1   );
μᵒᵘᵗ  = reshape(CuArray(-μ),  1,   n_μ, 1,  1,   1   );
_dϕ   = reshape(CuArray(dϕ), 1   ,1,   n_ϕ,1,   1   );
_μᴸ   = reshape(CuArray(μᴸ), 1   ,1,   1,  n_μᴸ,1   );
_dϕᴸ  = reshape(CuArray(dϕᴸ),1  , 1,   1,  1,   n_ϕᴸ);

cpu_μⁱⁿ   = reshape(μ,  n_μ, 1,   1,  1,   1   );
cpu_μᵒᵘᵗ  = reshape(-μ,  1,   n_μ, 1,  1,   1   );
cpu_dϕ    = reshape(dϕ, 1   ,1,   n_ϕ,1,   1   );
cpu_μᴸ    = reshape(μᴸ, 1   ,1,   1,  n_μᴸ,1   );
cpu_dϕᴸ   = reshape(dϕᴸ,1  , 1,   1,  1,   n_ϕᴸ);

@time CanopyOptics.leaf_dot_products.(μⁱⁿ, μᵒᵘᵗ, _dϕ, _μᴸ, _dϕᴸ);
@time CanopyOptics.leaf_dot_products.(cpu_μⁱⁿ, cpu_μᵒᵘᵗ, cpu_dϕ, cpu_μᴸ, cpu_dϕᴸ);

