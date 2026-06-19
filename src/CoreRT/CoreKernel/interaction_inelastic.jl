#=

This file contains RT interaction-related functions

=#

"""
Pre-allocated workspace for ScatteringInterface_11 interaction to avoid
repeated GPU memory allocation (~7 GB per call for typical Raman runs).
Create once before the layer loop and pass to `interaction!`.

When `staged=true`, uses CPU-backed output with per-pass GPU buffers:
- Only 3 GPU 4D output buffers at a time (vs 6), saving ~3.5 GB FP64
- Results staged through CPU between passes
"""
mutable struct InteractionWorkspace{A3, A4_mat, A4_src, A3_src}
    tmp_inv::A3         # nQuad × nQuad × nSpec
    tmpR⁻⁺::A3
    tmpR⁺⁻::A3
    tmpT⁻⁻::A3
    tmpT⁺⁺::A3
    tmpJ₀⁻::A3_src     # nQuad × 1 × nSpec
    tmpJ₀⁺::A3_src
    # GPU 4D output buffers (reused per pass: 2 mat + 1 src)
    gpu_ie_mat_A::A4_mat    # nQuad × nQuad × nSpec × nRaman
    gpu_ie_mat_B::A4_mat
    gpu_ie_src::A4_src      # nQuad × 1 × nSpec × nRaman
    # CPU-backed results for per-pass staging
    cpu_ieR⁻⁺::Array       # Pass 1 output
    cpu_ieT⁻⁻::Array
    cpu_ieJ₀⁻::Array
    cpu_ieR⁺⁻::Array       # Pass 2 output
    cpu_ieT⁺⁺::Array
    cpu_ieJ₀⁺::Array
    staged::Bool            # Whether to use per-pass CPU staging
end

function InteractionWorkspace(composite_layer, added_layer; staged::Bool=false)
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer
    @unpack ieR⁻⁺, ieR⁺⁻, ieT⁻⁻, ieT⁺⁺, ieJ₀⁺, ieJ₀⁻ = composite_layer
    # GPU 3D buffers (always needed)
    ws_tmp_inv = similar(T⁺⁺)
    ws_tmpR⁻⁺ = similar(R⁻⁺); ws_tmpR⁺⁻ = similar(R⁺⁻)
    ws_tmpT⁻⁻ = similar(T⁻⁻); ws_tmpT⁺⁺ = similar(T⁺⁺)
    ws_tmpJ₀⁻ = similar(J₀⁻); ws_tmpJ₀⁺ = similar(J₀⁺)
    # GPU 4D buffers: only 2 matrix + 1 source (reused per pass)
    ws_gpu_mat_A = similar(ieR⁻⁺)
    ws_gpu_mat_B = similar(ieR⁻⁺)
    ws_gpu_src   = similar(ieJ₀⁻)
    # CPU output arrays for staging
    FT = eltype(ieR⁻⁺)
    sz_mat = size(ieR⁻⁺)
    sz_src = size(ieJ₀⁻)
    InteractionWorkspace(
        ws_tmp_inv,
        ws_tmpR⁻⁺, ws_tmpR⁺⁻, ws_tmpT⁻⁻, ws_tmpT⁺⁺,
        ws_tmpJ₀⁻, ws_tmpJ₀⁺,
        ws_gpu_mat_A, ws_gpu_mat_B, ws_gpu_src,
        zeros(FT, sz_mat...), zeros(FT, sz_mat...), zeros(FT, sz_src...),  # Pass 1 CPU
        zeros(FT, sz_mat...), zeros(FT, sz_mat...), zeros(FT, sz_src...),  # Pass 2 CPU
        staged)
end

function reset!(ws::InteractionWorkspace)
    ws.tmp_inv .= 0
    ws.tmpR⁻⁺ .= 0; ws.tmpR⁺⁻ .= 0
    ws.tmpT⁻⁻ .= 0; ws.tmpT⁺⁺ .= 0
    ws.tmpJ₀⁻ .= 0; ws.tmpJ₀⁺ .= 0
    ws.gpu_ie_mat_A .= 0; ws.gpu_ie_mat_B .= 0; ws.gpu_ie_src .= 0
    return nothing
end

# No scattering in either the added layer or the composite layer
function interaction_helper!(RS_type,::ScatteringInterface_00, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    # If SFI, interact source function in no scattering
    # @show "interaction 00"
    if SFI
        composite_layer.ieJ₀⁺[:] = 0.0 #ieJ₀⁺
        composite_layer.ieJ₀⁻[:] = 0.0 #ieJ₀⁻

        composite_layer.J₀⁺ = added_layer.j₀⁺ .+ added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺
        composite_layer.J₀⁻ = composite_layer.J₀⁻ .+ composite_layer.T⁻⁻ ⊠ added_layer.j₀⁻
    end

    # Batched multiplication between added and composite
    composite_layer.T⁻⁻[:] = added_layer.t⁻⁻ ⊠ composite_layer.T⁻⁻
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺

    composite_layer.ieR⁻⁺[:] = 0.0
    composite_layer.ieR⁺⁻[:] = 0.0
    composite_layer.ieT⁻⁻[:] = 0.0
    composite_layer.ieT⁺⁺[:] = 0.0

end

# No scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer, added to bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper!(RS_type::RRS, ::ScatteringInterface_01, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    @unpack i_λ₁λ₀ = RS_type     
    
    # @show "interaction 01"
    if SFI   
        for n₁ in eachindex ieJ₁⁺[1,1,:,1]
            for Δn in eachindex ieJ₁⁺[1,1,1,:]
                n₀  = n₁ + i_λ₁λ₀[Δn]
                composite_layer.ieJ₀⁻[:,1,n₁,Δn] = 
                    composite_layer.T⁻⁻[:,:,n₁] * 
                    (added_layer.ier⁻⁺[:,:,n₁,Δn] * 
                    composite_layer.J₀⁺[:,1,n₀] + 
                    added_layer.ieJ₀⁻[:,1,n₁,Δn]) 
                composite_layer.ieJ₀⁺[:,1,n₁,Δn] = 
                    added_layer.ieJ₀⁺[:,1,n₁,Δn] + 
                    added_layer.iet⁺⁺[:,:,n₁,Δn] * 
                    composite_layer.J₀⁺[:,1,n₀]
            end 
        end

        composite_layer.J₀⁻ = composite_layer.J₀⁻ .+ composite_layer.T⁻⁻ ⊠ 
            (added_layer.r⁻⁺ ⊠ composite_layer.J₀⁺ .+ added_layer.j₀⁻) 
        composite_layer.J₀⁺ = added_layer.j₀⁺ .+ added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺ 
    end

    for n₁ in eachindex ieJ₁⁺[1,1,:,1]
        for Δn in eachindex ieJ₁⁺[1,1,1,:]
            n₀  = n₁ + i_λ₁λ₀[Δn]
            # Batched multiplication between added and composite
            composite_layer.ieR⁻⁺[:,:,n₁,Δn] = 
                composite_layer.T⁻⁻[:,:,n₁] * added_layer.ier⁻⁺[:,:,n₁,Δn] * 
                composite_layer.T⁺⁺[:,:,n₀]
    
            composite_layer.ieR⁺⁻[:,:,n₁,Δn] = added_layer.ier⁺⁻[:,:,n₁,Δn]
    
            composite_layer.ieT⁺⁺[:,:,n₁,Δn] = 
                added_layer.iet⁺⁺[:,:,n₁,Δn] * composite_layer.T⁺⁺[:,:,n₀]
    
            composite_layer.ieT⁻⁻[:,:,n₁,Δn] = 
                composite_layer.T⁻⁻[:,:,n₁] * added_layer.iet⁻⁻[:,:,n₁,Δn]    
        end
    end
    # Batched multiplication between added and composite
    composite_layer.R⁻⁺[:] = composite_layer.T⁻⁻ ⊠ added_layer.r⁻⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.R⁺⁻[:] = added_layer.r⁺⁻
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.T⁻⁻[:] = composite_layer.T⁻⁻ ⊠ added_layer.t⁻⁻    
end

function interaction_helper!(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
                        ::ScatteringInterface_01, SFI,
                        composite_layer::CompositeLayer{FT}, 
                        added_layer::AddedLayer{FT}, 
                        I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    @unpack i_λ₁λ₀_all = RS_type                             
    if SFI
        for Δn = 1:length(i_λ₁λ₀_all)
            n₁ = i_λ₁λ₀_all[Δn]
            n₀ = 1
            if (n₁>0)
                # TODO: replace all the n₀ indices with ref_r, ref_t, and ref_J
                composite_layer.ieJ₀⁻[:,1,n₁,n₀] = 
                    composite_layer.T⁻⁻[:,:,n₁] * 
                    (added_layer.ier⁻⁺[:,:,n₁,n₀] * 
                    composite_layer.J₀⁺[:,1,n₀] + 
                    added_layer.ieJ₀⁻[:,1,n₁,n₀]) 

                composite_layer.ieJ₀⁺[:,1,n₁,n₀] = 
                    added_layer.ieJ₀⁺[:,1,n₁,n₀] + 
                    added_layer.iet⁺⁺[:,:,n₁,n₀] * composite_layer.J₀⁺[:,1,n₀] 
            end
        end

        composite_layer.J₀⁻ = composite_layer.J₀⁻ .+ composite_layer.T⁻⁻ ⊠ 
            (added_layer.r⁻⁺ ⊠ composite_layer.J₀⁺ .+ added_layer.j₀⁻) 
        composite_layer.J₀⁺ = added_layer.j₀⁺ .+ added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺ 
    end

    for Δn = 1:length(i_λ₁λ₀_all)
        n₁ = i_λ₁λ₀_all[Δn]
        n₀ = 1
        # Batched multiplication between added and composite
        if (n₁>0)
            composite_layer.ieR⁻⁺[:,:,n₁,n₀] = 
                        composite_layer.T⁻⁻[:,:,n₁] * added_layer.ier⁻⁺[:,:,n₁,n₀] * 
                        composite_layer.T⁺⁺[:,:,n₀]

            composite_layer.ieR⁺⁻[:,:,n₁,n₀] = added_layer.ier⁺⁻[:,:,n₁,n₀]

            composite_layer.ieT⁺⁺[:,:,n₁,n₀] = 
                        added_layer.iet⁺⁺[:,:,n₁,n₀] * composite_layer.T⁺⁺[:,:,n₀]

            composite_layer.ieT⁻⁻[:,:,n₁,n₀] = 
                        composite_layer.T⁻⁻[:,:,n₁] * added_layer.iet⁻⁻[:,:,n₁,n₀] 
        end   
    end
    # Batched multiplication between added and composite
    composite_layer.R⁻⁺[:] = composite_layer.T⁻⁻ ⊠ added_layer.r⁻⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.R⁺⁻[:] = added_layer.r⁺⁻
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.T⁻⁻[:] = composite_layer.T⁻⁻ ⊠ added_layer.t⁻⁻    
end

# Scattering in inhomogeneous composite layer.
# no scattering in homogeneous layer which is 
# added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper!(RS_type::RRS, ::ScatteringInterface_10, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack i_λ₁λ₀ = RS_type 

    # @show "interaction 10"
    if SFI
        for n₁ in eachindex ieJ₁⁺[1,1,:,1]
            for Δn in eachindex ieJ₁⁺[1,1,1,:]
                n₀  = n₁ + i_λ₁λ₀[Δn]
                
                composite_layer.ieJ₀⁺[:,1,n₁,Δn] = 
                        added_layer.t⁺⁺[:,:,n₁] * 
                        (composite_layer.ieJ₀⁺[:,1,n₁,Δn] + 
                        composite_layer.ieR⁺⁻[:,:,n₁,Δn] * 
                        added_layer.j₀⁻[:,1,n₀])
                composite_layer.ieJ₀⁻[:,1,n₁,Δn] = 
                        composite_layer.ieJ₀⁻[:,1,n₁,Δn] + 
                        composite_layer.ieT⁻⁻[:,:,n₁,Δn] * 
                        added_layer.j₀⁻[:,1,n₀]
            end
        end
        composite_layer.J₀⁺ = added_layer.j₀⁺ .+ 
            added_layer.t⁺⁺ ⊠ (composite_layer.J₀⁺ .+ 
                composite_layer.R⁺⁻ ⊠ added_layer.j₀⁻)
        composite_layer.J₀⁻ = composite_layer.J₀⁻ .+ 
            composite_layer.T⁻⁻ ⊠ added_layer.j₀⁻
    end

    for n₁ in eachindex ieJ₁⁺[1,1,:,1]
        for Δn in eachindex ieJ₁⁺[1,1,1,:]
            n₀  = n₁ + i_λ₁λ₀[Δn]
            # Batched multiplication between added and composite
            composite_layer.ieT⁺⁺[:,:,n₁,Δn] = 
                    added_layer.t⁺⁺[:,:,n₁] * composite_layer.ieT⁺⁺[:,:,n₁,Δn]
            composite_layer.ieT⁻⁻[:,:,n₁,Δn] = 
                    composite_layer.ieT⁻⁻[:,:,n₁,Δn] * added_layer.t⁻⁻[:,:,n₀]
            composite_layer.ieR⁺⁻[:,:,n₁,Δn] = 
                    added_layer.t⁺⁺[:,:,n₁] * composite_layer.ieR⁺⁻[:,:,n₁,Δn] * 
                    added_layer.t⁻⁻[:,:,n₀]
        end
    end
    # Batched multiplication between added and composite
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.T⁻⁻[:] = composite_layer.T⁻⁻ ⊠ added_layer.t⁻⁻
    composite_layer.R⁺⁻[:] = added_layer.t⁺⁺ ⊠ composite_layer.R⁺⁻ ⊠ added_layer.t⁻⁻
end

function interaction_helper!(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
    ::ScatteringInterface_10, SFI,
    composite_layer::CompositeLayer{FT}, 
    added_layer::AddedLayer{FT}, 
    I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack i_λ₁λ₀_all = RS_type 
    if SFI
        for Δn = 1:length(i_λ₁λ₀_all)
            n₁ = i_λ₁λ₀_all[Δn]
            n₀ = 1
            if (n₁>0)
                composite_layer.ieJ₀⁺[:,1,n₁,n₀] = 
                    added_layer.t⁺⁺[:,:,n₁] * 
                    (composite_layer.ieJ₀⁺[:,1,n₁,n₀] + 
                    composite_layer.ieR⁺⁻[:,:,n₁,n₀] * 
                    added_layer.j₀⁻[:,1,n₀])

                composite_layer.ieJ₀⁻[:,1,n₁,n₀] = 
                            composite_layer.ieJ₀⁻[:,1,n₁,n₀] + 
                            composite_layer.ieT⁻⁻[:,:,n₁,n₀] * 
                            added_layer.j₀⁻[:,1,n₀]
            end
        end

        composite_layer.J₀⁺ = added_layer.j₀⁺ .+ 
            added_layer.t⁺⁺ ⊠ (composite_layer.J₀⁺ .+ 
            composite_layer.R⁺⁻ ⊠ added_layer.j₀⁻)
        composite_layer.J₀⁻ = composite_layer.J₀⁻ .+ 
            composite_layer.T⁻⁻ ⊠ added_layer.j₀⁻
    end

    for Δn = 1:length(i_λ₁λ₀_all)
        n₁ = i_λ₁λ₀_all[Δn]
        n₀ = 1
        # Batched multiplication between added and composite
        if (n₁>0)
            composite_layer.ieT⁺⁺[:,:,n₁,n₀] = 
                    added_layer.t⁺⁺[:,:,n₁] * composite_layer.ieT⁺⁺[:,:,n₁,n₀]
            composite_layer.ieT⁻⁻[:,:,n₁,n₀] = 
                composite_layer.ieT⁻⁻[:,:,n₁,n₀] * added_layer.t⁻⁻[:,:,n₀]
            composite_layer.ieR⁺⁻[:,:,n₁,n₀] = 
                added_layer.t⁺⁺[:,:,n₁] * composite_layer.ieR⁺⁻[:,:,n₁,n₀] * 
                added_layer.t⁻⁻[:,:,n₀]
        end
    end

    # Batched multiplication between added and composite
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.T⁻⁻[:] = composite_layer.T⁻⁻ ⊠ added_layer.t⁻⁻
    composite_layer.R⁺⁻[:] = added_layer.t⁺⁺ ⊠ composite_layer.R⁺⁻ ⊠ added_layer.t⁻⁻
end

# Scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer which is added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper!(RS_type::RRS, ::ScatteringInterface_11, SFI,
                                composite_layer::Union{CompositeLayer, CompositeLayerRS},
                                added_layer::Union{AddedLayer,AddedLayerRS},
                                I_static::AbstractArray{FT};
                                workspace::Union{InteractionWorkspace, Nothing}=nothing) where {FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    @unpack i_λ₁λ₀ = RS_type
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ = added_layer
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer
    @unpack ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺ = added_layer
    @unpack ieR⁻⁺, ieR⁺⁻, ieT⁺⁺, ieT⁻⁻, ieJ₀⁺, ieJ₀⁻ = composite_layer

    # Fused per-Δn matrix outputs (tmpieR⁻⁺/T⁻⁻ loop1, tmpieT⁺⁺/R⁺⁻ loop2) on GPU
    # when they fit static-@localmem/threads; else (CPU, or NquadN too large) the
    # batched_mul host loops below run. Source-vector outputs (tmpieJ₀±) stay
    # batched_mul. Fused path is not bit-identical → regen goldens for such configs.
    arch_ie = architecture(ier⁻⁺)
    use_fused_ie = (arch_ie isa Architectures.GPU) && _fused_interaction_fits(size(ier⁻⁺,1), eltype(ier⁻⁺))
    use_fused_ijsrc = (arch_ie isa Architectures.GPU) && _fused_ijsrc_fits(size(ier⁻⁺,1), eltype(ier⁻⁺))
    i_λ₁λ₀_dev = (use_fused_ie || use_fused_ijsrc) ? array_type(arch_ie)(i_λ₁λ₀) : i_λ₁λ₀

    # @show "interaction 11"
    staged = workspace !== nothing && workspace.staged
    # Used to store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    if workspace !== nothing
        reset!(workspace)
        tmp_inv = workspace.tmp_inv
        tmpJ₀⁻ = workspace.tmpJ₀⁻; tmpJ₀⁺ = workspace.tmpJ₀⁺
        tmpR⁻⁺ = workspace.tmpR⁻⁺; tmpR⁺⁻ = workspace.tmpR⁺⁻
        tmpT⁻⁻ = workspace.tmpT⁻⁻; tmpT⁺⁺ = workspace.tmpT⁺⁺
        if staged
            # Per-pass: use 3 shared GPU buffers, alias for Pass 1 outputs
            tmpieR⁻⁺ = workspace.gpu_ie_mat_A
            tmpieT⁻⁻ = workspace.gpu_ie_mat_B
            tmpieJ₀⁻ = workspace.gpu_ie_src
            # Pass 2 outputs will reuse same buffers after Pass 1 copies to CPU
            tmpieT⁺⁺ = workspace.gpu_ie_mat_A  # placeholder, reassigned after Pass 1
            tmpieR⁺⁻ = workspace.gpu_ie_mat_B
            tmpieJ₀⁺ = workspace.gpu_ie_src
        else
            tmpieJ₀⁻ = workspace.gpu_ie_src;  workspace.gpu_ie_src .= 0
            tmpieJ₀⁺ = similar(ieJ₀⁺); tmpieJ₀⁺ .= 0
            tmpieR⁻⁺ = workspace.gpu_ie_mat_A; workspace.gpu_ie_mat_A .= 0
            tmpieR⁺⁻ = workspace.gpu_ie_mat_B; workspace.gpu_ie_mat_B .= 0
            tmpieT⁻⁻ = similar(ieT⁻⁻); tmpieT⁻⁻ .= 0
            tmpieT⁺⁺ = similar(ieT⁺⁺); tmpieT⁺⁺ .= 0
        end
    else
        tmp_inv  = similar(t⁺⁺); tmp_inv.=0;
        tmpieJ₀⁻ = similar(ieJ₀⁻); tmpieJ₀⁻.=0;
        tmpieJ₀⁺ = similar(ieJ₀⁺); tmpieJ₀⁺.=0;
        tmpieR⁻⁺ = similar(ieR⁻⁺); tmpieR⁻⁺.=0;
        tmpieR⁺⁻ = similar(ieR⁺⁻); tmpieR⁺⁻.=0;
        tmpieT⁻⁻ = similar(ieT⁻⁻); tmpieT⁻⁻.=0;
        tmpieT⁺⁺ = similar(ieT⁺⁺); tmpieT⁺⁺.=0;
        tmpJ₀⁻   = similar(J₀⁻); tmpJ₀⁻.=0;
        tmpJ₀⁺   = similar(J₀⁺); tmpJ₀⁺.=0;
        tmpR⁻⁺   = similar(R⁻⁺); tmpR⁻⁺.=0;
        tmpR⁺⁻   = similar(R⁺⁻); tmpR⁺⁻.=0;
        tmpT⁻⁻   = similar(T⁻⁻); tmpT⁻⁻.=0;
        tmpT⁺⁺   = similar(T⁺⁺); tmpT⁺⁺.=0;
        staged = false
    end
    # Compute and store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    @timeit "interaction inv1" batch_inv!(tmp_inv, I_static .- r⁻⁺ ⊠ R⁺⁻) #Suniti
    # Temporary arrays:
    # T₁₂(I-R₀₁R₂₁)⁻¹
    T01_inv = T⁻⁻ ⊠ tmp_inv;
    if SFI
        if use_fused_ijsrc
            apply_fused_interaction_Jm!(tmpieJ₀⁻, T01_inv, tmp_inv, ier⁻⁺, r⁻⁺, R⁺⁻, ieR⁺⁻, ieT⁻⁻, ieJ₀⁻, J₀⁺, ieJ₀⁺, added_layer.ieJ₀⁻, added_layer.j₀⁻, i_λ₁λ₀_dev, arch_ie)
        else
        for Δn = 1:size(ieJ₀⁺,4)
            n₀, n₁ = get_n₀_n₁(ieJ₀⁺,i_λ₁λ₀[Δn])
            @inbounds @views tmpieJ₀⁻[:,:,n₁,Δn] .=
                    ieJ₀⁻[:,:,n₁,Δn] +
                    T01_inv[:,:,n₁] ⊠ 
                    (ier⁻⁺[:,:,n₁,Δn] ⊠ J₀⁺[:,:,n₀] +
                    r⁻⁺[:,:,n₁] ⊠ ieJ₀⁺[:,:,n₁,Δn] +
                    added_layer.ieJ₀⁻[:,:,n₁,Δn]) +
                    # ── Suniti, please re-derive. The (n₀, n₁, Δn) index
                    #    triple here was the site of a long-standing
                    #    sign/index mismatch in the inelastic interaction
                    #    update. Verify that the second `T01_inv ⊠ (…)`
                    #    block below indexes the recipient wavelength `n₁`
                    #    and donor wavelength `n₀` consistently with the
                    #    Hovenier composition rule used in `interaction!`
                    #    for the elastic part.
                    (T01_inv[:,:,n₁] ⊠
                    (ier⁻⁺[:,:,n₁,Δn] ⊠ R⁺⁻[:,:,n₀] + 
                    r⁻⁺[:,:,n₁] ⊠ ieR⁺⁻[:,:,n₁,Δn]) +
                    ieT⁻⁻[:,:,n₁,Δn]) ⊠
                    tmp_inv[:,:,n₀] ⊠ 
                    (added_layer.j₀⁻[:,:,n₀] + r⁻⁺[:,:,n₀] ⊠ J₀⁺[:,:,n₀]);
        end
        end  # if use_fused_ijsrc (tmpieJ₀⁻)
        #J₀₂⁻ = J₀₁⁻ + T₀₁(1-R₂₁R₀₁)⁻¹(R₂₁J₁₀⁺+J₁₂⁻)
        tmpJ₀⁻ .= J₀⁻ .+ T01_inv ⊠ (r⁻⁺ ⊠ J₀⁺ .+ added_layer.j₀⁻) 
    end 
    if use_fused_ie
        apply_fused_interaction_Rmp!(tmpieR⁻⁺, T01_inv, tmp_inv, r⁻⁺, R⁺⁻, T⁺⁺, ier⁻⁺, ieR⁺⁻, ieT⁺⁺, ieT⁻⁻, ieR⁻⁺, i_λ₁λ₀_dev, arch_ie)
        apply_fused_interaction_Tmm!(tmpieT⁻⁻, T01_inv, tmp_inv, r⁻⁺, R⁺⁻, t⁻⁻, ier⁻⁺, ieR⁺⁻, ieT⁻⁻, iet⁻⁻, i_λ₁λ₀_dev, arch_ie)
    else
    for Δn = 1:size(ier⁻⁺,4)
        n₀, n₁ = get_n₀_n₁(ier⁻⁺,i_λ₁λ₀[Δn])
        @inbounds @views tmpieR⁻⁺[:,:,n₁,Δn] .= ieR⁻⁺[:,:,n₁,Δn] +
            T01_inv[:,:,n₁] ⊠   
            (ier⁻⁺[:,:,n₁,Δn] ⊠ T⁺⁺[:,:,n₀] + r⁻⁺[:,:,n₁] ⊠ ieT⁺⁺[:,:,n₁,Δn]) +    
            (T01_inv[:,:,n₁] ⊠ 
            (ier⁻⁺[:,:,n₁,Δn] ⊠ R⁺⁻[:,:,n₀] + r⁻⁺[:,:,n₁] ⊠ ieR⁺⁻[:,:,n₁,Δn]) + 
            ieT⁻⁻[:,:,n₁,Δn]) ⊠ 
            tmp_inv[:,:,n₀] ⊠ r⁻⁺[:,:,n₀] ⊠ T⁺⁺[:,:,n₀]
        @inbounds @views tmpieT⁻⁻[:,:,n₁,Δn] .=
            T01_inv[:,:,n₁] ⊠ iet⁻⁻[:,:,n₁,Δn] +  
            (T01_inv[:,:,n₁] ⊠ 
            (ier⁻⁺[:,:,n₁,Δn] ⊠ R⁺⁻[:,:,n₀] + r⁻⁺[:,:,n₁] ⊠ ieR⁺⁻[:,:,n₁,Δn]) +
            ieT⁻⁻[:,:,n₁,Δn]) ⊠ 
            tmp_inv[:,:,n₀] ⊠ t⁻⁻[:,:,n₀]
    end
    end  # if use_fused_ie (loop 1: tmpieR⁻⁺, tmpieT⁻⁻)

    # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀
    tmpR⁻⁺ .= R⁻⁺ .+ T01_inv ⊠ r⁻⁺ ⊠ T⁺⁺ #Suniti
    # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
    tmpT⁻⁻ .= T01_inv ⊠ t⁻⁻ #Suniti

    # --- Stage Pass 1 results to CPU, prepare Pass 2 GPU buffers ---
    if staged
        synchronize_if_gpu()
        copyto!(workspace.cpu_ieR⁻⁺, tmpieR⁻⁺)  # gpu_ie_mat_A → CPU
        copyto!(workspace.cpu_ieT⁻⁻, tmpieT⁻⁻)  # gpu_ie_mat_B → CPU
        copyto!(workspace.cpu_ieJ₀⁻, tmpieJ₀⁻)  # gpu_ie_src → CPU
        workspace.gpu_ie_mat_A .= 0
        workspace.gpu_ie_mat_B .= 0
        workspace.gpu_ie_src .= 0
        tmpieT⁺⁺ = workspace.gpu_ie_mat_A
        tmpieR⁺⁻ = workspace.gpu_ie_mat_B
        tmpieJ₀⁺ = workspace.gpu_ie_src
    end

    # Repeating for mirror-reflected directions

    # Compute and store `(I - r⁻⁺ * R⁺⁻)⁻¹`
    @timeit "interaction inv2" batch_inv!(tmp_inv, I_static .- R⁺⁻ ⊠ r⁻⁺) #Suniti
    # T₂₁(I-R₀₁R₂₁)⁻¹
    T21_inv = t⁺⁺ ⊠ tmp_inv
    if SFI
        if use_fused_ijsrc
            apply_fused_interaction_Jp!(tmpieJ₀⁺, T21_inv, tmp_inv, ier⁻⁺, r⁻⁺, R⁺⁻, ieR⁺⁻, iet⁺⁺, ieJ₀⁺, J₀⁺, added_layer.ieJ₀⁺, added_layer.ieJ₀⁻, added_layer.j₀⁻, i_λ₁λ₀_dev, arch_ie)
        else
        for Δn = 1:size(ieJ₀⁺,4)
            n₀, n₁ = get_n₀_n₁(ieJ₀⁺,i_λ₁λ₀[Δn])
            @inbounds @views tmpieJ₀⁺[:,:,n₁,Δn] .=
                            added_layer.ieJ₀⁺[:,:,n₁,Δn] +
                            T21_inv[:,:,n₁] ⊠ 
                            (ieJ₀⁺[:,:,n₁,Δn] + 
                            ieR⁺⁻[:,:,n₁,Δn] ⊠ added_layer.j₀⁻[:,:,n₀] +
                            R⁺⁻[:,:,n₁] ⊠ added_layer.ieJ₀⁻[:,:,n₁,Δn]) +
                            (T21_inv[:,:,n₁] ⊠ 
                            (ieR⁺⁻[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀] + 
                            R⁺⁻[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn]) +
                            iet⁺⁺[:,:,n₁,Δn]) ⊠ 
                            tmp_inv[:,:,n₀] ⊠ (J₀⁺[:,:,n₀] + 
                            R⁺⁻[:,:,n₀] ⊠ added_layer.j₀⁻[:,:,n₀])
        end
        end  # if use_fused_ijsrc (tmpieJ₀⁺)
        # J₂₀⁺ = J₂₁⁺ + T₂₁(I-R₀₁R₂₁)⁻¹(J₁₀ + R₀₁J₁₂⁻ )
        tmpJ₀⁺ = added_layer.j₀⁺ .+ 
            T21_inv ⊠ (J₀⁺ + R⁺⁻ ⊠ added_layer.j₀⁻)
    end 
    if use_fused_ie
        apply_fused_interaction_Tpp!(tmpieT⁺⁺, T21_inv, tmp_inv, r⁻⁺, R⁺⁻, T⁺⁺, ier⁻⁺, ieR⁺⁻, ieT⁺⁺, iet⁺⁺, i_λ₁λ₀_dev, arch_ie)
        apply_fused_interaction_Rpm2!(tmpieR⁺⁻, T21_inv, tmp_inv, r⁻⁺, R⁺⁻, t⁻⁻, ier⁻⁺, ieR⁺⁻, iet⁻⁻, iet⁺⁺, ier⁺⁻, i_λ₁λ₀_dev, arch_ie)
    else
    for Δn = 1:size(ieJ₀⁺,4)
        n₀, n₁ = get_n₀_n₁(ieJ₀⁺,i_λ₁λ₀[Δn])

        @inbounds @views tmpieT⁺⁺[:,:,n₁,Δn] .=
                    T21_inv[:,:,n₁] ⊠ ieT⁺⁺[:,:,n₁,Δn] +
                    (T21_inv[:,:,n₁] ⊠ (ieR⁺⁻[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀] + 
                    R⁺⁻[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn]) +
                    iet⁺⁺[:,:,n₁,Δn]) ⊠ tmp_inv[:,:,n₀] ⊠ T⁺⁺[:,:,n₀] #Suniti: Eq 12 of Raman paper draft

        @inbounds @views tmpieR⁺⁻[:,:,n₁,Δn] .=
                    ier⁺⁻[:,:,n₁,Δn] +
                    T21_inv[:,:,n₁] ⊠ 
                    (ieR⁺⁻[:,:,n₁,Δn] ⊠ t⁻⁻[:,:,n₀] +
                    R⁺⁻[:,:,n₁] ⊠ iet⁻⁻[:,:,n₁,Δn]) +
                    (T21_inv[:,:,n₁] ⊠ 
                    (ieR⁺⁻[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀] + 
                    R⁺⁻[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn]) + 
                    iet⁺⁺[:,:,n₁,Δn]) ⊠
                    tmp_inv[:,:,n₀] ⊠ R⁺⁻[:,:,n₀] ⊠ t⁻⁻[:,:,n₀]
    end
    end  # if use_fused_ie (loop 2: tmpieT⁺⁺, tmpieR⁺⁻)

    # T₂₀ = T₂₁(I-R₀₁R₂₁)⁻¹T₁₀
    tmpT⁺⁺ .= T21_inv  ⊠ T⁺⁺ #Suniti
    # R₀₂ = R₁₂ + T₂₁(1-R₀₁R₂₁)⁻¹R₀₁T₁₂
    tmpR⁺⁻ .= r⁺⁻ .+ T21_inv ⊠ R⁺⁻ ⊠ t⁻⁻ #Suniti

    composite_layer.J₀⁻ .= tmpJ₀⁻
    composite_layer.R⁻⁺ .= tmpR⁻⁺
    composite_layer.T⁻⁻ .= tmpT⁻⁻

    composite_layer.J₀⁺ .= tmpJ₀⁺
    composite_layer.T⁺⁺ .= tmpT⁺⁺
    composite_layer.R⁺⁻ .= tmpR⁺⁻

    if staged
        # Stage Pass 2 results to CPU
        synchronize_if_gpu()
        copyto!(workspace.cpu_ieT⁺⁺, tmpieT⁺⁺)  # gpu_ie_mat_A → CPU
        copyto!(workspace.cpu_ieR⁺⁻, tmpieR⁺⁻)  # gpu_ie_mat_B → CPU
        copyto!(workspace.cpu_ieJ₀⁺, tmpieJ₀⁺)  # gpu_ie_src → CPU
        # Copy ALL staged results from CPU → GPU composite_layer fields
        copyto!(composite_layer.ieJ₀⁻, workspace.cpu_ieJ₀⁻)
        copyto!(composite_layer.ieR⁻⁺, workspace.cpu_ieR⁻⁺)
        copyto!(composite_layer.ieT⁻⁻, workspace.cpu_ieT⁻⁻)
        copyto!(composite_layer.ieJ₀⁺, workspace.cpu_ieJ₀⁺)
        copyto!(composite_layer.ieT⁺⁺, workspace.cpu_ieT⁺⁺)
        copyto!(composite_layer.ieR⁺⁻, workspace.cpu_ieR⁺⁻)
    else
        composite_layer.ieJ₀⁻ .= tmpieJ₀⁻
        composite_layer.ieJ₀⁺ .= tmpieJ₀⁺
        composite_layer.ieT⁻⁻ .= tmpieT⁻⁻
        composite_layer.ieR⁻⁺ .= tmpieR⁻⁺
        composite_layer.ieT⁺⁺ .= tmpieT⁺⁺
        composite_layer.ieR⁺⁻ .= tmpieR⁺⁻
    end
end

function interaction_helper!(RS_type::Union{VS_0to1_plus, VS_1to0_plus},
    ::ScatteringInterface_11, SFI,
    composite_layer::Union{CompositeLayer, CompositeLayerRS},
    added_layer::Union{AddedLayer,AddedLayerRS},
    I_static::AbstractArray{FT};
    workspace::Union{InteractionWorkspace, Nothing}=nothing) where {FT<:Union{AbstractFloat, ForwardDiff.Dual}}

    @unpack i_λ₁λ₀_all = RS_type
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺ = added_layer
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻, ieR⁻⁺, ieR⁺⁻, ieT⁺⁺, ieT⁻⁻, ieJ₀⁺, ieJ₀⁻ = composite_layer
    # Used to store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    if workspace !== nothing
        reset!(workspace)
        tmp_inv  = workspace.tmp_inv
        tmpieJ₀⁻ = workspace.tmpieJ₀⁻; tmpieJ₀⁺ = workspace.tmpieJ₀⁺
        tmpieR⁻⁺ = workspace.tmpieR⁻⁺; tmpieR⁺⁻ = workspace.tmpieR⁺⁻
        tmpieT⁻⁻ = workspace.tmpieT⁻⁻; tmpieT⁺⁺ = workspace.tmpieT⁺⁺
        tmpJ₀⁻ = workspace.tmpJ₀⁻; tmpJ₀⁺ = workspace.tmpJ₀⁺
        tmpR⁻⁺ = workspace.tmpR⁻⁺; tmpR⁺⁻ = workspace.tmpR⁺⁻
        tmpT⁻⁻ = workspace.tmpT⁻⁻; tmpT⁺⁺ = workspace.tmpT⁺⁺
    else
        tmp_inv = similar(t⁺⁺); tmp_inv.=0;
        tmpieJ₀⁻ = similar(ieJ₀⁻); tmpieJ₀⁻.=0;
        tmpieJ₀⁺ = similar(ieJ₀⁺); tmpieJ₀⁺.=0;
        tmpieR⁻⁺ = similar(ieR⁻⁺); tmpieR⁻⁺.=0;
        tmpieR⁺⁻ = similar(ieR⁺⁻); tmpieR⁺⁻.=0;
        tmpieT⁻⁻ = similar(ieT⁻⁻); tmpieT⁻⁻.=0;
        tmpieT⁺⁺ = similar(ieT⁺⁺); tmpieT⁺⁺.=0;
        tmpJ₀⁻   = similar(J₀⁻); tmpJ₀⁻.=0;
        tmpJ₀⁺   = similar(J₀⁺); tmpJ₀⁺.=0;
        tmpR⁻⁺   = similar(R⁻⁺); tmpR⁻⁺.=0;
        tmpR⁺⁻   = similar(R⁺⁻); tmpR⁺⁻.=0;
        tmpT⁻⁻   = similar(T⁻⁻); tmpT⁻⁻.=0;
        tmpT⁺⁺   = similar(T⁺⁺); tmpT⁺⁺.=0;
    end
    # Compute and store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    @timeit "interaction inv1" batch_inv!(tmp_inv, I_static .- r⁻⁺ ⊠ R⁺⁻) #Suniti
    # Temporary arrays:
    # T₁₂(I-R₀₁R₂₁)⁻¹
    T01_inv = T⁻⁻ ⊠ tmp_inv;

    if SFI
        for Δn = 1:length(i_λ₁λ₀_all) # in eachindex ieJ₁⁺[1,1,:,1]
            n₁ = i_λ₁λ₀_all[Δn]
            n₀ = 1
            if (n₁>0)
                @inbounds @views tmpieJ₀⁻[:,1,n₁,n₀] = ieJ₀⁻[:,1,n₁,n₀] + 
                                        T01_inv[:,:,n₁] *
                                        (ier⁻⁺[:,:,n₁,n₀] * J₀⁺[:,1,n₀] + 
                                        r⁻⁺[:,:,n₁] * ieJ₀⁺[:,1,n₁,n₀] +
                                        added_layer.ieJ₀⁻[:,1,n₁,n₀]) +
                                        (T01_inv[:,:,n₁] * 
                                        (ier⁻⁺[:,:,n₁,n₀] * R⁺⁻[:,:,n₀] + 
                                        r⁻⁺[:,:,n₁] * ieR⁺⁻[:,:,n₁,n₀]) +
                                        ieT⁻⁻[:,:,n₁,n₀]) *
                                        tmp_inv[:,:,n₀] * 
                                        (added_layer.j₀⁻[:,1,n₀] + r⁻⁺[:,:,n₀] * J₀⁺[:,1,n₀]) #Suniti: Eq 17 of Raman paper draft
            end
        end
        #J₀₂⁻ = J₀₁⁻ + T₀₁(1-R₂₁R₀₁)⁻¹(R₂₁J₁₀⁺+J₁₂⁻)
        tmpJ₀⁻ .= J₀⁻ .+ T01_inv ⊠ (r⁻⁺ ⊠ J₀⁺ .+ added_layer.j₀⁻) 
    end 

    for Δn = 1:length(i_λ₁λ₀_all) # in eachindex ieJ₁⁺[1,1,:,1]
        n₁ = i_λ₁λ₀_all[Δn]
        n₀ = 1
        if (n₁>0)
            @inbounds @views tmpieR⁻⁺[:,:,n₁,n₀] = ieR⁻⁺[:,:,n₁,n₀] +
                    T01_inv[:,:,n₁] * 
                    (ier⁻⁺[:,:,n₁,n₀] * T⁺⁺[:,:,n₀] + r⁻⁺[:,:,n₁] * ieT⁺⁺[:,:,n₁,n₀]) +    
                    (T01_inv[:,:,n₁] * 
                    (ier⁻⁺[:,:,n₁,n₀] * R⁺⁻[:,:,n₀] + r⁻⁺[:,:,n₁] * ieR⁺⁻[:,:,n₁,n₀]) + 
                    ieT⁻⁻[:,:,n₁,n₀]) * 
                    tmp_inv[:,:,n₀] * r⁻⁺[:,:,n₀] * T⁺⁺[:,:,n₀] #Suniti: Eq 14 of Raman paper draft

            @inbounds @views tmpieT⁻⁻[:,:,n₁,n₀] = T01_inv[:,:,n₁] * iet⁻⁻[:,:,n₁,n₀] +  
                    (T01_inv[:,:,n₁] * 
                    (ier⁻⁺[:,:,n₁,n₀] * R⁺⁻[:,:,n₀] + r⁻⁺[:,:,n₁] * ieR⁺⁻[:,:,n₁,n₀]) +
                    ieT⁻⁻[:,:,n₁,n₀]) * 
                    tmp_inv[:,:,n₀] * t⁻⁻[:,:,n₀] #Suniti: Eq 13 of Raman paper draft
        end
    end

    # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀ 
    tmpR⁻⁺ .= R⁻⁺ .+ T01_inv ⊠ r⁻⁺ ⊠ T⁺⁺ #Suniti
    # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
    tmpT⁻⁻ .= T01_inv ⊠ t⁻⁻ #Suniti

    # Repeating for mirror-reflected directions

    # Compute and store `(I - r⁻⁺ * R⁺⁻)⁻¹`
    @timeit "interaction inv2" batch_inv!(tmp_inv, I_static .- R⁺⁻ ⊠ r⁻⁺) #Suniti
    # T₂₁(I-R₀₁R₂₁)⁻¹
    T21_inv = t⁺⁺ ⊠ tmp_inv
    if SFI
        for Δn = 1:length(i_λ₁λ₀_all) #Δn in eachindex ieJ₁⁺[1,1,:,1]
            n₁ = i_λ₁λ₀_all[Δn]
            n₀ = 1
            if (n₁>0)
                tmpieJ₀⁺[:,1,n₁,n₀] = added_layer.ieJ₀⁺[:,1,n₁,n₀] + 
                        T21_inv[:,:,n₁] * 
                        (ieJ₀⁺[:,1,n₁,n₀] + 
                        ieR⁺⁻[:,:,n₁,n₀] * added_layer.j₀⁻[:,1,n₀] +
                        R⁺⁻[:,:,n₁] * added_layer.ieJ₀⁻[:,1,n₁,n₀]) +
                        (T21_inv[:,:,n₁] * 
                        (ieR⁺⁻[:,:,n₁,n₀] * r⁻⁺[:,:,n₀] + 
                        R⁺⁻[:,:,n₁] * ier⁻⁺[:,:,n₁,n₀]) +
                        iet⁺⁺[:,:,n₁,n₀]) * 
                        tmp_inv[:,:,n₀] * (J₀⁺[:,1,n₀] + 
                        R⁺⁻[:,:,n₀] * added_layer.j₀⁻[:,1,n₀])
            end
        end
    # J₂₀⁺ = J₂₁⁺ + T₂₁(I-R₀₁R₂₁)⁻¹(J₁₀ + R₀₁J₁₂⁻ )
    tmpJ₀⁺ = added_layer.j₀⁺ .+ 
                T21_inv ⊠ (J₀⁺ .+ R⁺⁻ ⊠ added_layer.j₀⁻)
    end 
    for Δn = 1:length(i_λ₁λ₀_all) #Δn in eachindex ieJ₁⁺[1,1,:,1]
        n₁ = i_λ₁λ₀_all[Δn]
        n₀ = 1
        if (n₁>0)
            tmpieT⁺⁺[:,:,n₁,n₀] = T21_inv[:,:,n₁] * ieT⁺⁺[:,:,n₁,n₀] +
                    (T21_inv[:,:,n₁] * (ieR⁺⁻[:,:,n₁,n₀] * r⁻⁺[:,:,n₀] + 
                    R⁺⁻[:,:,n₁] * ier⁻⁺[:,:,n₁,n₀]) +
                    iet⁺⁺[:,:,n₁,n₀]) * tmp_inv[:,:,n₀] * T⁺⁺[:,:,n₀] #Suniti: Eq 12 of Raman paper draft

            tmpieR⁺⁻[:,:,n₁,n₀] = ier⁺⁻[:,:,n₁,n₀] + 
                    T21_inv[:,:,n₁] *
                    (ieR⁺⁻[:,:,n₁,n₀] * t⁻⁻[:,:,n₀] + 
                    R⁺⁻[:,:,n₁] * iet⁻⁻[:,:,n₁,n₀]) +
                    (T21_inv[:,:,n₁] * 
                    (ieR⁺⁻[:,:,n₁,n₀] * r⁻⁺[:,:,n₀] + 
                    R⁺⁻[:,:,n₁] * ier⁻⁺[:,:,n₁,n₀]) + 
                    iet⁺⁺[:,:,n₁,n₀]) *
                    tmp_inv[:,:,n₀] * R⁺⁻[:,:,n₀] * t⁻⁻[:,:,n₀] #Suniti: Eq 15 of Raman paper draft
        end
    end
    # T₂₀ = T₂₁(I-R₀₁R₂₁)⁻¹T₁₀
    tmpT⁺⁺ .= T21_inv  ⊠ T⁺⁺ #Suniti
    # R₀₂ = R₁₂ + T₂₁(1-R₀₁R₂₁)⁻¹R₀₁T₁₂
    tmpR⁺⁻ .= r⁺⁻ .+ T21_inv ⊠ R⁺⁻ ⊠ t⁻⁻ #Suniti

    composite_layer.J₀⁻ .= tmpJ₀⁻
    composite_layer.R⁻⁺ .= tmpR⁻⁺
    composite_layer.T⁻⁻ .= tmpT⁻⁻

    composite_layer.J₀⁺ .= tmpJ₀⁺
    composite_layer.T⁺⁺ .= tmpT⁺⁺
    composite_layer.R⁺⁻ .= tmpR⁺⁻
    
    composite_layer.ieJ₀⁻ .= tmpieJ₀⁻
    composite_layer.ieJ₀⁺ .= tmpieJ₀⁺

    composite_layer.ieT⁻⁻ .= tmpieT⁻⁻
    composite_layer.ieR⁻⁺ .= tmpieR⁻⁺
    composite_layer.ieT⁺⁺ .= tmpieT⁺⁺
    composite_layer.ieR⁺⁻ .= tmpieR⁺⁻

end

"Compute interaction between composite and added layers"
function interaction!(RS_type::Union{RRS, VS_0to1_plus, VS_1to0_plus}, scattering_interface::AbstractScatteringInterface, SFI,
                        composite_layer::Union{CompositeLayer,CompositeLayerRS},
                        added_layer::Union{AddedLayer,AddedLayerRS},
                        I_static::AbstractArray{FT};
                        workspace::Union{InteractionWorkspace, Nothing}=nothing) where {FT<:Union{AbstractFloat, ForwardDiff.Dual}}

    interaction_helper!(RS_type, scattering_interface, SFI, composite_layer, added_layer, I_static; workspace=workspace)
    
    #M1 = composite_layer.R⁻⁺[1,1,1]
    #M2 = composite_layer.R⁺⁻[1,1,1]
    #M3 = composite_layer.T⁺⁺[1,1,1]
    #M4 = composite_layer.T⁻⁻[1,1,1]
    #M5 = composite_layer.J₀⁺[1,1,1]
    #M6 = composite_layer.J₀⁻[1,1,1]

    #@show M1, M2, M3, M4, M5, M6
    synchronize_if_gpu()
    
end

function interaction!(RS_type::noRS, scattering_interface::AbstractScatteringInterface, SFI,
    composite_layer::Union{CompositeLayer,CompositeLayerRS},
    added_layer::Union{AddedLayer,AddedLayerRS},
    I_static::AbstractArray{FT};
    workspace=nothing) where {FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    
    interaction_helper!(scattering_interface, SFI, composite_layer, added_layer, I_static)
    synchronize_if_gpu() # @timeit "interaction" boundary — ensures device work completes before timing stops

end