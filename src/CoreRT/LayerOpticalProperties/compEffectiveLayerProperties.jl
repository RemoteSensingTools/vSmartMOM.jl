function constructCoreOpticalProperties(RS_type::AbstractRamanType, iBand, m, model)
    (; τ_rayl, τ_aer, τ_abs, aerosol_optics, greek_rayleigh, greek_cabannes) = model
    @assert all(iBand .≤ length(τ_rayl)) "iBand exceeded number of bands"

    arr_type = CoreRT.array_type(model)

    pol_type = CoreRT.polarization_type(model)

    # Quadrature points:
    μ = collect(model.quad_points.qp_μ)
    N = length(model.quad_points.qp_μN)
    # Number of Aerosols:
    nAero = size(τ_aer[iBand[1]],1)
    nZ    = size(τ_rayl[1],2)
    # Rayleigh Z matrix per band — pick greek_rayleigh for noRS (pure-elastic,
    # includes rotational Raman in the effective depol) or greek_cabannes for
    # any Raman-aware RS_type (rotational Raman is handled separately via the
    # RRS interaction, so the elastic path should use the lower-depol Cabannes
    # phase matrix). Mismatch here produces a ~1% bias on Stokes I and ~3% on Q
    # because the Cabannes depol (~0.007) vs Rayleigh depol (~0.028) shifts the
    # polarization-sensitive greek coefficients (β, δ) by ~3%.
    # See sanghavi src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl:31-42.
    raman_active = !(RS_type isa noRS) && !(RS_type isa noRS_plus)
    band_layer_props    = Vector{Vector{CoreScatteringOpticalProperties}}()
    band_fScattRayleigh = Vector{Vector}()
    for iB in iBand
        gr_source = raman_active ? greek_cabannes : greek_rayleigh
        gr = gr_source isa AbstractVector ? gr_source[iB] : gr_source
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, μ, gr, m,
                                                        arr_type = arr_type)
        rayl = [CoreScatteringOpticalProperties(arr_type(τ_rayl[iB][:,i]), RS_type.ϖ_Cabannes[iB],
                Rayl𝐙⁺⁺, Rayl𝐙⁻⁺) for i=1:nZ]
        
        combo = rayl

        for i=1:nAero
            AerZ⁺⁺, AerZ⁻⁺ = Scattering.compute_Z_moments(
                                pol_type, μ, 
                                aerosol_optics[iB][i].greek_coefs, 
                                m, arr_type=arr_type)
            aer = createAero.(τ_aer[iB][i,:], 
                              [aerosol_optics[iB][i]], 
                              [AerZ⁺⁺], [AerZ⁻⁺])
            combo = combo .+ aer
        end

        fScattRayleigh = [Array(rayl[i].τ ./ combo[i].τ) for i=1:nZ]
        push!(band_layer_props,
              combo .+ [CoreAbsorptionOpticalProperties(arr_type((τ_abs[iB][:,i]))) for i=1:nZ])
        push!(band_fScattRayleigh, fScattRayleigh)
    end

    layer_opt = [prod([band_layer_props[i][iz] for i=1:length(iBand)]) for iz=1:nZ]
    fscat_opt = [[band_fScattRayleigh[i][iz] for i=1:length(iBand)] for iz=1:nZ]
    return layer_opt, fscat_opt
end

function createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺)
    (; fᵗ, ω̃) = aerosol_optics
    τ_mod = (1-fᵗ * ω̃ ) * τAer;
    ϖ_mod = (1-fᵗ) * ω̃/(1-fᵗ * ω̃)
    CoreScatteringOpticalProperties(τ_mod, ϖ_mod,AerZ⁺⁺, AerZ⁻⁺)
end

# Extract scattering definitions and integrated absorptions for the source function!
function extractEffectiveProps(
            lods::Array,#{CoreScatteringOpticalProperties{FT},1}
            quad_points::QuadPoints{FT}
            ) where FT

    nSpec = length(lods[1].τ)
    nZ    = length(lods)
    scattering_interface = ScatteringInterface_00()
    scattering_interfaces_all = AbstractScatteringInterface[]
    τ_sum_all = similar(lods[1].τ, (nSpec, nZ+1))
    τ_sum_all[:,1] .= 0
    for iz = 1:nZ
        scatter = maximum(lods[iz].τ .* lods[iz].ϖ) > 2eps(FT)
        scattering_interface = get_scattering_interface(scattering_interface, scatter, iz)
        push!(scattering_interfaces_all, scattering_interface)
        @views τ_sum_all[:,iz+1] = τ_sum_all[:,iz] + getG_atSun(lods[iz], quad_points) * lods[iz].τ 
    end
    return scattering_interfaces_all, τ_sum_all
end

function getG_atSun(lod::CoreScatteringOpticalProperties,quad_points::QuadPoints{FT}) where FT
    FT(1)
end

function getG_atSun(lod::CoreDirectionalScatteringOpticalProperties,quad_points::QuadPoints{FT}) where FT
    (; iμ₀) = quad_points
    gfct = collect(lod.G)[iμ₀]
    return gfct
end


function expandOpticalProperties(in::CoreScatteringOpticalProperties, arr_type)
    (; τ, ϖ, Z⁺⁺, Z⁻⁺) = in 
    @assert length(τ) == length(ϖ) "τ and ϖ sizes need to match"
    if size(Z⁺⁺,3) == 1
        Z⁺⁺ = _repeat(Z⁺⁺,1,1,length(τ))
        Z⁻⁺ = _repeat(Z⁻⁺,1,1,length(τ))
        return CoreScatteringOpticalProperties(arr_type(τ), arr_type(ϖ), arr_type(Z⁺⁺), arr_type(Z⁻⁺)) 
    else
        @assert size(Z⁺⁺,3) ==  length(τ) "Z and τ dimensions need to match "
        CoreScatteringOpticalProperties(arr_type(τ), arr_type(ϖ), arr_type(Z⁺⁺), arr_type(Z⁻⁺)) 
    end
end

function expandOpticalProperties(in::CoreDirectionalScatteringOpticalProperties, arr_type)
    (; τ, ϖ, Z⁺⁺, Z⁻⁺, G) = in 
    @assert length(τ) == length(ϖ) "τ and ϖ sizes need to match"
    if size(Z⁺⁺,3) == 1
        Z⁺⁺ = _repeat(Z⁺⁺,1,1,length(τ))
        Z⁻⁺ = _repeat(Z⁻⁺,1,1,length(τ))
        return CoreDirectionalScatteringOpticalProperties(arr_type(τ), arr_type(ϖ), arr_type(Z⁺⁺), arr_type(Z⁻⁺), arr_type(G)) 
    else
        @assert size(Z⁺⁺,3) ==  length(τ) "Z and τ dimensions need to match "
        CoreDirectionalScatteringOpticalProperties(arr_type(τ), arr_type(ϖ), arr_type(Z⁺⁺), arr_type(Z⁻⁺), arr_type(G)) 
    end
end

function expandBandScalars(RS_type, x)
    out = zeros(eltype(x[1]), sum([length(RS_type.bandSpecLim[iB]) for iB in RS_type.iBand]))
    for iB in RS_type.iBand
        out[RS_type.bandSpecLim[iB]] .= expandScalar(x[iB],length(RS_type.bandSpecLim[iB]))
    end
    return out
end

expandScalar(x,n) = x.*ones(n);