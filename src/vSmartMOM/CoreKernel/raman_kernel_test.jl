nij = 24
nn1 = 1000
nn0 = 100

iet⁺⁺ = rand(nij,nij,nn1,nn0);
ier⁻⁺ = rand(nij,nij,nn1,nn0);
using KernelAbstractions

ϖ_λ₀λ₁ = rand(nn1,nn0);
dτ₀= rand(1)[1];
dτ₁= rand(1)[1];
dτ_λ= rand(nn1);
Z⁻⁺_λ₀λ₁= rand(nij,nij);
Z⁺⁺_λ₀λ₁= rand(nn1,nn0);
qp_μN= rand(nij);
wct2 = rand(nij);
ϖ_λ = rand(nn1);

@kernel function get_elem_rt!(ier⁻⁺, iet⁺⁺, ϖ_λ, ϖ_λ₀λ₁, dτ₀, dτ₁, dτ_λ, Z⁻⁺_λ₀λ₁, Z⁺⁺_λ₀λ₁, qp_μN, wct2)
    i, j, n₁, n₀ = @index(Global, NTuple) 
    # let n₁ cover the full range of wavelengths, while n₀ only includes wavelengths at intervals 
    # that contribute significantly enough to inelastic scattering, so that n₀≪n₁ 
    if (wct2[j]>1.e-8) 
        # dτ₀, dτ₁ are the purely scattering (elastic+inelastic) molecular elemental 
        # optical thicknesses at wavelengths λ₀ and λ₁
        # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        ier⁻⁺[i,j,n₁,n₀] = ϖ_λ₀λ₁[n₁,n₀] * (dτ₀/dτ₁) * Z⁻⁺_λ₀λ₁[i,j] * (qp_μN[j]*dτ₁ / (qp_μN[i]*dτ₀ + qp_μN[j]*dτ₁)) * (1 - exp(-((dτ_λ[n₁] / qp_μN[i]) + (dτ_λ[n₀] / qp_μN[j])))) * (wct2[j]) 
                    
        if (qp_μN[i] == qp_μN[j])
            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
            if i == j       
                if abs(dτ_λ[n₀]-dτ_λ[n₁])>1.e-6
                    iet⁺⁺[i,j,n₁,n₀] = ((exp(-dτ_λ[n₀] / qp_μN[i]) - exp(-dτ_λ[n₁] / qp_μN[i]))/(dτ_λ[n₁]-dτ_λ[n₀])) * ϖ_λ₀λ₁[n₁,n₀] * dτ₀ * Z⁺⁺_λ₀λ₁[i,i] * wct2[i]
                else    
                    iet⁺⁺[i,j,n₁,n₀] = ϖ_λ₀λ₁[n₁,n₀] * dτ₀ * Z⁺⁺_λ₀λ₁[i,i] * wct2[i] * exp(-dτ_λ[n₀] / qp_μN[j])/ qp_μN[j]
                end
            else
                iet⁺⁺[i,j,n₁,n₀] = 0.0
            end
        else  
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)
            iet⁺⁺[i,j,n₁,n₀] = ϖ_λ₀λ₁[n₁,n₀] * (dτ₀/dτ₁) * Z⁺⁺_λ₀λ₁[i,j] * (qp_μN[j]*dτ₁ / (qp_μN[i]*dτ₀ - qp_μN[j]*dτ₁)) * (exp(-dτ_λ[n₁] / qp_μN[i]) - exp(-dτ_λ[n₀] / qp_μN[j])) * wct2[j]
        end
    else
        ier⁻⁺[i,j,n₁,n₀] = 0.0
        if i==j
            iet⁺⁺[i,j,n₁,n₀] = 0.0
        else
            iet⁺⁺[i,j,n₁,n₀] = 0.0
        end
    end
end

function get_elem_rt!(ier⁻⁺, iet⁺⁺, ϖ_λ, ϖ_λ₀λ₁, dτ₀, dτ₁, dτ_λ, Z⁻⁺_λ₀λ₁, Z⁺⁺_λ₀λ₁, qp_μN, wct2)
    d1,d2,d3,d4 = size(ier⁻⁺)
    for i in 1:d1, j in 1:d2, n₁ in 1:d3, n₀ in 1:d4
        #i, j, n₁, n₀ = @index(Global, NTuple) 
        # let n₁ cover the full range of wavelengths, while n₀ only includes wavelengths at intervals 
        # that contribute significantly enough to inelastic scattering, so that n₀≪n₁ 
        if (wct2[j]>1.e-8) 
            # dτ₀, dτ₁ are the purely scattering (elastic+inelastic) molecular elemental 
            # optical thicknesses at wavelengths λ₀ and λ₁
            # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
            ier⁻⁺[i,j,n₁,n₀] = ϖ_λ₀λ₁[n₁,n₀] * (dτ₀/dτ₁) * Z⁻⁺_λ₀λ₁[i,j] * (qp_μN[j]*dτ₁ / (qp_μN[i]*dτ₀ + qp_μN[j]*dτ₁)) * (1 - exp(-((dτ_λ[n₁] / qp_μN[i]) + (dτ_λ[n₀] / qp_μN[j])))) * (wct2[j]) 
                        
            if (qp_μN[i] == qp_μN[j])
                # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
                if i == j       
                    if abs(dτ_λ[n₀]-dτ_λ[n₁])>1.e-6
                        iet⁺⁺[i,j,n₁,n₀] = ((exp(-dτ_λ[n₀] / qp_μN[i]) - exp(-dτ_λ[n₁] / qp_μN[i]))/(dτ_λ[n₁]-dτ_λ[n₀])) * ϖ_λ₀λ₁[n₁,n₀] * dτ₀ * Z⁺⁺_λ₀λ₁[i,i] * wct2[i]
                    else    
                        iet⁺⁺[i,j,n₁,n₀] = ϖ_λ₀λ₁[n₁,n₀] * dτ₀ * Z⁺⁺_λ₀λ₁[i,i] * wct2[i] * exp(-dτ_λ[n₀] / qp_μN[j])/ qp_μN[j]
                    end
                else
                    iet⁺⁺[i,j,n₁,n₀] = 0.0
                end
            else  
                # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
                # (𝑖 ≠ 𝑗)
                iet⁺⁺[i,j,n₁,n₀] = ϖ_λ₀λ₁[n₁,n₀] * (dτ₀/dτ₁) * Z⁺⁺_λ₀λ₁[i,j] * (qp_μN[j]*dτ₁ / (qp_μN[i]*dτ₀ - qp_μN[j]*dτ₁)) * (exp(-dτ_λ[n₁] / qp_μN[i]) - exp(-dτ_λ[n₀] / qp_μN[j])) * wct2[j]
            end
        else
            ier⁻⁺[i,j,n₁,n₀] = 0.0
            if i==j
                iet⁺⁺[i,j,n₁,n₀] = 0.0
            else
                iet⁺⁺[i,j,n₁,n₀] = 0.0
            end
        end
    end
end
# Test w/o kernel:
get_elem_rt!(ier⁻⁺, iet⁺⁺, ϖ_λ, ϖ_λ₀λ₁,dτ₀, dτ₁, dτ_λ, Z⁻⁺_λ₀λ₁, Z⁺⁺_λ₀λ₁, qp_μN, wct2);

# Test CPU kernel version:
device = CPU()
kernel! = get_elem_rt!(device)
event = kernel!(ier⁻⁺, iet⁺⁺, ϖ_λ, ϖ_λ₀λ₁,dτ₀, dτ₁, dτ_λ, Z⁻⁺_λ₀λ₁, Z⁺⁺_λ₀λ₁, qp_μN, wct2, ndrange=size(ier⁻⁺)); 
wait(device, event)