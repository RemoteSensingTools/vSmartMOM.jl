
using Test
using NNlib
import NNlib.batched_mul
using CUDA

"Batched matrix multiply (overwrite NNlib definition)"
function batched_mul(A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    CUBLAS.gemm_strided_batched('N', 'N', A, B)
end

nSpec = 2000;
nQuad = 10;
nRaman = 134;
aType = CuArray

ieJ₁⁺ = aType(rand(nQuad, 1,nSpec,nRaman));
ieJ₁⁻ = aType(rand(nQuad, 1,nSpec,nRaman));

ieJ₀⁻ = aType(rand(nQuad, 1,nSpec,nRaman));
ieJ₀⁺ = aType(rand(nQuad, 1,nSpec,nRaman));
ieJ₀⁻_ = deepcopy(ieJ₀⁻);
ieJ₀⁺_ = deepcopy(ieJ₀⁺);
ier⁻⁺ = aType(rand(nQuad, nQuad,nSpec,nRaman));
iet⁺⁺ = aType(rand(nQuad, nQuad,nSpec,nRaman));
iet⁻⁻ = aType(rand(nQuad, nQuad,nSpec,nRaman));
r⁻⁺   = aType(rand(nQuad, nQuad,nSpec));
tt⁺⁺_gp_refl = aType(rand(nQuad, nQuad,nSpec));
gp_refl = aType(rand(nQuad, nQuad,nSpec));
J₁⁻ = aType(rand(nQuad, 1,nSpec));
J₀⁺ = aType(rand(nQuad, 1,nSpec));
i_λ₁λ₀ = collect(1:nRaman);


 # in eachindex ieJ₁⁺[1,1,:,1]
 @views function compLoop!(ieJ₀⁺,ieJ₀⁻,ieJ₁⁺, ieJ₁⁻, ier⁻⁺, iet⁺⁺, iet⁻⁻, r⁻⁺, tt⁺⁺_gp_refl,J₁⁻, J₀⁺, i_λ₁λ₀ ) 
    for n₁ = 1:size(ieJ₁⁺,3)
        for Δn = 1:size(ieJ₁⁺,4) #in eachindex ieJ₁⁺[1,1,1,:]

            n₀  = n₁ + i_λ₁λ₀[Δn]
            if 1 ≤ n₀ ≤ size(ieJ₁⁺,3)
                # J⁺₂₀(λ) = J⁺₂₁(λ) + T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹[J⁺₁₀(λ) + R⁺⁻₀₁(λ)J⁻₁₂(λ)] (see Eqs.16 in Raman paper draft)
                @inbounds ieJ₀⁺[:,1,n₁,Δn] = ieJ₁⁺[:,1,n₁,Δn] + 
                        (tt⁺⁺_gp_refl[:,:,n₁] * 
                        (ieJ₀⁺[:,1,n₁,Δn] + r⁻⁺[:,:,n₁] * ieJ₁⁻[:,1,n₁,Δn] + ier⁻⁺[:,:,n₁,Δn] * J₁⁻[:,1,n₀] + 
                        (r⁻⁺[:,:,n₁] * ier⁻⁺[:,:,n₁,Δn] + ier⁻⁺[:,:,n₁,Δn] * r⁻⁺[:,:,n₀]) * 
                        gp_refl[:,:,n₀] * (J₀⁺[:,1,n₀] + r⁻⁺[:,:,n₀] * J₁⁻[:,1,n₀]))) + 
                        iet⁺⁺[:,:,n₁,Δn] * gp_refl[:,:,n₀] * 
                        (J₀⁺[:,1,n₀] + r⁻⁺[:,:,n₀] * J₁⁻[:,1,n₀])
                # J⁻₀₂(λ) = J⁻₀₁(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹[J⁻₁₂(λ) + R⁻⁺₂₁(λ)J⁺₁₀(λ)] (see Eqs.17 in Raman paper draft)
                @inbounds ieJ₀⁻[:,1,n₁,Δn] = ieJ₀⁻[:,1,n₁,Δn] + 
                        (tt⁺⁺_gp_refl[:,:,n₁] * 
                        (ieJ₁⁻[:,1,n₁,Δn] + ier⁻⁺[:,:,n₁,Δn] * J₀⁺[:,1,n₀] + r⁻⁺[:,:,n₁] * ieJ₀⁺[:,1,n₁,Δn] + 
                        (ier⁻⁺[:,:,n₁,Δn] * r⁻⁺[:,:,n₀] + r⁻⁺[:,:,n₁] * ier⁻⁺[:,:,n₁,Δn]) * 
                        gp_refl[:,:,n₀] * (J₁⁻[:,1,n₀] + r⁻⁺[:,:,n₀] * J₀⁺[:,1,n₀]))) +
                        iet⁻⁻[:,:,n₁,Δn] * gp_refl[:,:,n₀] * (J₁⁻[:,1,n₀] + 
                        r⁻⁺[:,:,n₀]*J₀⁺[:,1,n₀])
            end
        end
    end 
    return nothing
end

function batched_mul(A::AbstractArray{FT,3}, B::AbstractArray{FT,4}) where {FT}
    return reshape(NNlib.batched_mul(A,reshape(B,(size(B,1),size(B,2), size(B,3)*size(B,4)))), size(B))
end

function compBatched!(ieJ₀⁺,ieJ₀⁻,ieJ₁⁺, ieJ₁⁻, ier⁻⁺, iet⁺⁺, iet⁻⁻, r⁻⁺, tt⁺⁺_gp_refl,J₁⁻, J₀⁺, i_λ₁λ₀ )
    # Version 2, batched:
    tmp1 = gp_refl ⊠  (J₀⁺ + r⁻⁺ ⊠ J₁⁻)
    tmp2 = gp_refl ⊠  (J₁⁻ + r⁻⁺ ⊠ J₀⁺)
    for Δn = 1:size(ieJ₁⁺,4)
        n₁_ = 1:size(ieJ₁⁺,3);
        n₀_  = n₁_ .+ i_λ₁λ₀[Δn];
        # Find valid indices:
        sub = findall(1 .≤ n₀_ .≤ size(ieJ₁⁺,3));
        n₁ = view(n₁_,sub);
        n₀ = view(n₀_,sub);

        #@show length(n₁), length(n₀), length(n₁_), length(n₀_)
        @inbounds ieJ₀⁺[:,:,n₁,Δn] = ieJ₁⁺[:,:,n₁,Δn] + 
                        (tt⁺⁺_gp_refl[:,:,n₁] ⊠ 
                        (ieJ₀⁺[:,:,n₁,Δn] + r⁻⁺[:,:,n₁] ⊠ ieJ₁⁻[:,:,n₁,Δn] + ier⁻⁺[:,:,n₁,Δn] ⊠ J₁⁻[:,:,n₀] + 
                        (r⁻⁺[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn] + ier⁻⁺[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀]) ⊠ 
                        tmp1[:,:,n₀])) + 
                        iet⁺⁺[:,:,n₁,Δn] ⊠ tmp1[:,:,n₀]
        @inbounds ieJ₀⁻[:,:,n₁,Δn] = ieJ₀⁻[:,:,n₁,Δn] + 
                        (tt⁺⁺_gp_refl[:,:,n₁] ⊠ 
                        (ieJ₁⁻[:,:,n₁,Δn] + ier⁻⁺[:,:,n₁,Δn] ⊠ J₀⁺[:,:,n₀] + r⁻⁺[:,:,n₁] ⊠ ieJ₀⁺[:,:,n₁,Δn] + 
                        (ier⁻⁺[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀] + r⁻⁺[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn]) ⊠ 
                        tmp2[:,:,n₀])) +
                        iet⁻⁻[:,:,n₁,Δn] ⊠ tmp2[:,:,n₀]
    end
    return nothing
end


@time    compLoop!(ieJ₀⁺,ieJ₀⁻,ieJ₁⁺, ieJ₁⁻, ier⁻⁺, iet⁺⁺, iet⁻⁻, r⁻⁺, tt⁺⁺_gp_refl,J₁⁻, J₀⁺, i_λ₁λ₀ )
@time compBatched!(ieJ₀⁺_,ieJ₀⁻_,ieJ₁⁺, ieJ₁⁻, ier⁻⁺, iet⁺⁺, iet⁻⁻, r⁻⁺, tt⁺⁺_gp_refl,J₁⁻, J₀⁺, i_λ₁λ₀ )

@test ieJ₀⁺ ≈ ieJ₀⁺_
@test ieJ₀⁻ ≈ ieJ₀⁻_

function gpu_MM(A, B, C)
    for i = 1:size(B,4)
        @inbounds C[:,:,:,i] = batched_mul(A,view(B,:,:,:,i)); 
    end
    return nothing
end
