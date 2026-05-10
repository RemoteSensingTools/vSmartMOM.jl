ier⁻⁺ = rand(20,20,1000,100);
iet⁻⁻ = rand(20,20,1000,100);
gp_refl = rand(20,20,1000);
ieJ₀⁺ = rand(20, 1, 1000,100);
ieJ₁⁻ = rand(20, 1, 1000,100);
r⁻⁺   = rand(20,20,1000);
J₁⁻= rand(20, 1, 1000);

function NNlib.batched_mul(A::Array{FT,3}, B::Vector{FT}) where {FT}
    return NNlib.batched_mul(A,reshape(B,(size(B,1),1)))
end

@kernel function get_doubling_ie_rtSFI!(RS_type::RRS, 
    r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻,
    ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺, ieJ₀⁺, ieJ₀⁻)

-, -, n₁, Δn = @index(Global, NTuple)
@unpack fscattRayl, ϖ_λ₁λ₀, i_λ₁λ₀, i_ref = RS_type 
# n₁ covers the full range of wavelengths, while n₀ = n₁+Δn only includes wavelengths at intervals 
# that contribute significantly enough to inelastic scattering, so that n₀≪n₁ 
n₀  = n₁ + i_RRS[Δn]

# J⁻₀₂(λ) = J⁻₀₁(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹[J⁻₁₂(λ) + R⁻⁺₂₁(λ)J⁺₁₀(λ)] (see Eqs.17 in Raman paper draft)
ieJ₀⁻[:,1,n₁,Δn] = ieJ₀⁻[:,1,n₁,Δn] .+ (tt⁺⁺_gp_refl[:,:,n₁] ⊠ 
    (ieJ₁⁻[:,1,n₁,Δn] .+
    ier⁻⁺[:,:,n₁,Δn] ⊠ J₀⁺[:,1,n₀] .+ 
    r⁻⁺[:,:,n₁] ⊠ ieJ₀⁺[:,1,n₁,Δn] .+ 
    (ier⁻⁺[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀] .+ r⁻⁺[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn]) ⊠ 
    gp_refl[:,:,n₀] ⊠ (J₁⁻[:,1,n₀] .+ r⁻⁺[:,:,n₀] ⊠ J₀⁺[:,1,n₀]))) .+
    iet⁻⁻[:,:,n₁,Δn] ⊠ gp_refl[:,:,n₀] ⊠ 
    (J₁⁻[:,1,n₀] .+ r⁻⁺[:,:,n₀] ⊠ J₀⁺[:,1,n₀])

# J⁺₂₀(λ) = J⁺₂₁(λ) + T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹[J⁺₁₀(λ) + R⁺⁻₀₁(λ)J⁻₁₂(λ)] (see Eqs.16 in Raman paper draft)
ieJ₀⁺[:,1,n₁,Δn] = ieJ₁⁺[:,1,n₁,Δn] .+ 
    (tt⁺⁺_gp_refl[:,:,n₁] ⊠ (ieJ₀⁺[:,1,n₁,Δn] .+ 
     r⁻⁺[:,:,n₁] ⊠ ieJ₁⁻[:,1,n₁,Δn] .+ ier⁻⁺[:,:,n₁,Δn] ⊠ J₁⁻[:,1,n₀] .+ 
     (r⁻⁺[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn] .+ ier⁻⁺[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀]) ⊠ 
     gp_refl[:,:,n₀] ⊠ (J₀⁺[:,1,n₀] .+ r⁻⁺[:,:,n₀] ⊠ J₁⁻[:,1,n₀]))) .+ 
     iet⁺⁺[:,:,n₁,Δn] ⊠ gp_refl[:,:,n₀] ⊠ 
     (J₀⁺[:,1,n₀] .+ r⁻⁺[:,:,n₀] ⊠ J₁⁻[:,1,n₀])

if (wct2[j]>1.e-8) 

# dτ₀, dτ₁ are the purely scattering (elastic+inelastic) molecular elemental 
# optical thicknesses at wavelengths λ₀ and λ₁
# 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
ier⁻⁺[i,j,n₁,Δn] = fscattRayl * ϖ_λ₁λ₀[i_ϖ] * Z⁻⁺_λ₁λ₀[i,j] * 
(qp_μN[j] / (qp_μN[i] + qp_μN[j])) * 
(1 - exp(-((dτ_λ[n₁] / qp_μN[i]) + (dτ_λ[n₀] / qp_μN[j])))) * wct2[j] 

if (qp_μN[i] == qp_μN[j])
# @show i,j
# 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
if i == j       
if abs(dτ_λ[n₀]-dτ_λ[n₁])>1.e-6
iet⁺⁺[i,j,n₁,Δn] = 
ϖ_λ₁λ₀[i_ϖ] * fscattRayl * dτ * Z⁺⁺_λ₁λ₀[i,i] * wct2[i] *
((exp(-dτ_λ[n₀] / qp_μN[i]) - exp(-dτ_λ[n₁] / qp_μN[i]))/(dτ_λ[n₁]-dτ_λ[n₀])) 

else    
iet⁺⁺[i,j,n₁,Δn] = 
ϖ_λ₁λ₀[i_ϖ] * fscattRayl * dτ * Z⁺⁺_λ₁λ₀[i,i] * wct2[i] *
exp(-dτ_λ[n₀] / qp_μN[j])/ qp_μN[j]
end
else
iet⁺⁺[i,j,n₁,Δn] = 0.0
end
else
#@show  qp_μN[i], qp_μN[j]  
# 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
# (𝑖 ≠ 𝑗)
iet⁺⁺[i,j,n₁,Δn] = 
ϖ_λ₁λ₀[i_ϖ] * fscattRayl * Z⁺⁺_λ₁λ₀[i,j] * 
(qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct2[j] * 
(exp(-dτ_λ[n₁] / qp_μN[i]) - exp(-dτ_λ[n₀] / qp_μN[j]))
end
else
ier⁻⁺[i,j,n₁,Δn] = 0.0
if i==j
iet⁺⁺[i,j,n₁,Δn] = 0.0
else
iet⁺⁺[i,j,n₁,Δn] = 0.0
end
end
end