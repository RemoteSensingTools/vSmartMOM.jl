#=

Pretty-printing for vSmartMOM core types (Oceananigans-inspired tree style)

=#

# ── Tree-drawing constants ────────────────────────────────────────────────
const _TREE_MID   = "├── "
const _TREE_END   = "└── "
const _TREE_PIPE  = "│   "
const _TREE_SPACE = "    "

# ── Surface one-line summaries ────────────────────────────────────────────

function Base.show(io::IO, s::LambertianSurfaceScalar)
    print(io, "LambertianSurfaceScalar(α=$(s.albedo))")
end

function Base.show(io::IO, s::LambertianSurfaceSpectrum)
    n = length(s.albedo)
    print(io, "LambertianSurfaceSpectrum($n points)")
end

function Base.show(io::IO, s::LambertianSurfaceLegendre)
    n = length(s.legendre_coeff)
    print(io, "LambertianSurfaceLegendre(degree=$(n-1))")
end

function Base.show(io::IO, s::LambertianSurfaceSpline)
    n = length(s.wlGrid)
    print(io, "LambertianSurfaceSpline($n knots)")
end

function Base.show(io::IO, s::rpvSurfaceScalar)
    print(io, "rpvSurface(ρ₀=$(s.ρ₀), k=$(s.k), Θ=$(s.Θ), ρ_c=$(s.ρ_c))")
end

function Base.show(io::IO, s::RossLiSurfaceScalar)
    print(io, "RossLi(fvol=$(s.fvol), fgeo=$(s.fgeo), fiso=$(s.fiso))")
end

function Base.show(io::IO, s::CoxMunkSurface)
    wc = s.include_whitecaps ? "whitecaps=on" : "whitecaps=off"
    print(io, "CoxMunkSurface(U=$(s.wind_speed) m/s, $wc)")
end

function Base.show(io::IO, s::CanopySurface)
    soil_name = typeof(s.soil).name.name
    print(io, "CanopySurface(LAI=$(s.LAI), $(s.n_layers) layer(s), soil=$soil_name)")
end

# ── Layer one-line summaries ──────────────────────────────────────────────

function Base.show(io::IO, l::CompositeLayer{FT}) where FT
    n1, n2, ns = size(l.R⁻⁺)
    print(io, "CompositeLayer{$FT}($(n1)×$(n2), $ns spec)")
end

function Base.show(io::IO, l::AddedLayer{FT}) where FT
    n1, n2, ns = size(l.r⁻⁺)
    print(io, "AddedLayer{$FT}($(n1)×$(n2), $ns spec)")
end

# ── Atmosphere ───────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", a::Atmosphere{FT}) where FT
    prof = a.profile
    nLevels = length(prof.p_full)
    nLayers = nLevels
    nBands = length(a.spec_bands)

    println(io, "Atmosphere{$FT}")

    # profile
    p_top = round(minimum(prof.p_half), sigdigits=4)
    p_bot = round(maximum(prof.p_half), sigdigits=4)
    println(io, _TREE_MID, "profile: $nLayers layers, p ∈ [$p_top, $p_bot] hPa")
    println(io, _TREE_PIPE, _TREE_MID, "T range: [$(round(minimum(prof.T), digits=1)), $(round(maximum(prof.T), digits=1))] K")
    # VMR gases
    gas_keys = collect(keys(prof.vmr))
    if !isempty(gas_keys)
        println(io, _TREE_PIPE, _TREE_END, "gases: ", join(sort(gas_keys), ", "))
    else
        println(io, _TREE_PIPE, _TREE_END, "gases: none")
    end

    # spectral bands
    band_strs = [begin
        b = a.spec_bands[i]
        "[$(round(b[1], digits=1))–$(round(b[end], digits=1))]"
    end for i in 1:min(nBands, 4)]
    band_summary = join(band_strs, " ")
    if nBands > 4
        band_summary *= " …"
    end
    print(io, _TREE_END, "spectral_bands: $nBands band(s) $band_summary cm⁻¹")
end

function Base.show(io::IO, a::Atmosphere{FT}) where FT
    nLayers = length(a.profile.p_full)
    nBands = length(a.spec_bands)
    print(io, "Atmosphere{$FT}($nLayers layers, $nBands band(s))")
end

# ── AtmosphericProfile ──────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", p::AtmosphericProfile{FT}) where FT
    nLevels = length(p.p_full)
    p_top = round(minimum(p.p_half), sigdigits=4)
    p_bot = round(maximum(p.p_half), sigdigits=4)

    println(io, "AtmosphericProfile{$FT}")
    println(io, _TREE_MID, "layers: $nLevels")
    println(io, _TREE_MID, "p range: [$p_top, $p_bot] hPa")
    println(io, _TREE_MID, "T range: [$(round(minimum(p.T), digits=1)), $(round(maximum(p.T), digits=1))] K")
    gas_keys = collect(keys(p.vmr))
    if !isempty(gas_keys)
        print(io, _TREE_END, "gases: ", join(sort(gas_keys), ", "))
    else
        print(io, _TREE_END, "gases: none")
    end
end

function Base.show(io::IO, p::AtmosphericProfile{FT}) where FT
    nLevels = length(p.p_full)
    print(io, "AtmosphericProfile{$FT}($nLevels layers)")
end

# ── Optics ───────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", o::Optics{FT}) where FT
    nBands = length(o.τ_abs)
    nAer = length(o.aerosols.aerosol_optics) > 0 ? length(o.aerosols.aerosol_optics[1]) : 0

    println(io, "Optics{$FT}")

    # Rayleigh
    nCab = length(o.rayleigh.ϖ_Cabannes)
    ϖ_range = if nCab == 1
        "$(round(o.rayleigh.ϖ_Cabannes[1], digits=4))"
    else
        "$(round(minimum(o.rayleigh.ϖ_Cabannes), digits=4))–$(round(maximum(o.rayleigh.ϖ_Cabannes), digits=4))"
    end
    println(io, _TREE_MID, "rayleigh: $nCab band(s), ϖ_Cabannes=$ϖ_range")

    # Aerosols
    if nAer > 0
        println(io, _TREE_MID, "aerosols: $nAer type(s)")
        for iB in 1:min(length(o.aerosols.aerosol_optics), 1)
            aer_ops = o.aerosols.aerosol_optics[iB]
            for (ia, ao) in enumerate(aer_ops)
                ω_str = ao.ω̃ isa Number ? "$(round(ao.ω̃, digits=4))" : "$(round(minimum(ao.ω̃), digits=4))–$(round(maximum(ao.ω̃), digits=4))"
                fᵗ_str = ao.fᵗ isa Number ? "$(round(ao.fᵗ, digits=4))" : "$(round(minimum(ao.fᵗ), digits=4))–$(round(maximum(ao.fᵗ), digits=4))"
                l_max = length(ao.greek_coefs.β) - 1
                connector = (ia < length(aer_ops)) ? _TREE_MID : _TREE_END
                println(io, _TREE_PIPE, connector, "aerosol #$ia: ω̃=$ω_str, fᵗ=$fᵗ_str, l_max=$l_max")
            end
        end
    else
        println(io, _TREE_MID, "aerosols: none")
    end

    # τ_abs
    if nBands > 0 && length(o.τ_abs) > 0
        τ_shape = size(o.τ_abs[1])
        println(io, _TREE_MID, "τ_abs: $nBands band(s), $(τ_shape[1])×$(τ_shape[2]) per band")
    end

    # τ_rayl
    if nBands > 0 && length(o.τ_rayl) > 0
        τ_shape = size(o.τ_rayl[1])
        print(io, _TREE_END, "τ_rayl: $nBands band(s), $(τ_shape[1])×$(τ_shape[2]) per band")
    end
end

function Base.show(io::IO, o::Optics{FT}) where FT
    nBands = length(o.τ_abs)
    nAer = length(o.aerosols.aerosol_optics) > 0 ? length(o.aerosols.aerosol_optics[1]) : 0
    print(io, "Optics{$FT}($nBands band(s), $nAer aerosol(s))")
end

# ── RTModel ───────────────────────────────────────────────────────────────

function Base.summary(io::IO, m::RTModel)
    FT = float_type(m)
    nBands = length(m.atmosphere.spec_bands)
    nLevels = length(m.atmosphere.profile.p_full)
    nStreams = m.quad_points.Nstreams
    print(io, "RTModel{$(typeof(m.architecture)), $FT}($nStreams weighted streams, $nBands band(s), $nLevels levels)")
end

function Base.show(io::IO, ::MIME"text/plain", m::RTModel)
    FT = float_type(m)
    nBands = length(m.atmosphere.spec_bands)
    nAer = n_aerosols(m)
    nLevels = length(m.atmosphere.profile.p_full)
    prof = m.atmosphere.profile

    println(io, "RTModel{$(typeof(m.architecture)), $FT}")

    # architecture
    println(io, _TREE_MID, "architecture: ", m.architecture)

    # solver
    pol_name = typeof(m.solver.polarization_type).name.name
    println(io, _TREE_MID, "solver: $pol_name, Nstreams=$(m.quad_points.Nstreams), Nquad=$(m.quad_points.Nquad), max_m=$(m.solver.max_m), l_trunc=$(m.solver.l_trunc)")

    # geometry
    vza_str = length(m.geometry.vza) <= 4 ? string(m.geometry.vza) : "[$(length(m.geometry.vza)) angles]"
    println(io, _TREE_MID, "geometry: SZA=$(m.geometry.sza)°, VZA=$(vza_str)°")

    # quad_points
    println(io, _TREE_MID, "quad_points: $(m.quad_points.Nstreams) weighted streams (Nquad=$(m.quad_points.Nquad) inc. SZA/VZA output nodes), μ₀=$(round(m.quad_points.μ₀, digits=4))")

    # atmosphere
    band_strs = [begin
        b = m.atmosphere.spec_bands[i]
        "[$(round(b[1], digits=1))–$(round(b[end], digits=1))]"
    end for i in 1:min(nBands, 4)]
    band_summary = join(band_strs, " ")
    if nBands > 4
        band_summary *= " …"
    end
    p_top = round(minimum(prof.p_half), sigdigits=4)
    p_bot = round(maximum(prof.p_half), sigdigits=4)
    println(io, _TREE_MID, "atmosphere: $nLevels layers, $nBands band(s)")
    println(io, _TREE_PIPE, _TREE_MID, "p range: [$p_top, $p_bot] hPa")
    println(io, _TREE_PIPE, _TREE_MID, "T range: [$(round(minimum(prof.T), digits=1)), $(round(maximum(prof.T), digits=1))] K")
    gas_keys = collect(keys(prof.vmr))
    if !isempty(gas_keys)
        println(io, _TREE_PIPE, _TREE_MID, "gases: ", join(sort(gas_keys), ", "))
    end
    println(io, _TREE_PIPE, _TREE_END, "bands: $band_summary cm⁻¹")

    # optics
    println(io, _TREE_MID, "optics: $nAer aerosol(s)")
    if nBands > 0 && length(m.optics.τ_abs) > 0
        τ_shape = size(m.optics.τ_abs[1])
        println(io, _TREE_PIPE, _TREE_MID, "τ_abs: $nBands band(s), $(τ_shape[1])×$(τ_shape[2]) per band")
    end
    # aerosol details
    if nAer > 0 && length(m.optics.aerosols.aerosol_optics) > 0
        for iB in 1:min(length(m.optics.aerosols.aerosol_optics), 1)  # show band 1 only
            aer_ops = m.optics.aerosols.aerosol_optics[iB]
            for (ia, ao) in enumerate(aer_ops)
                ω_str = ao.ω̃ isa Number ? "$(round(ao.ω̃, digits=4))" : "$(round(minimum(ao.ω̃), digits=4))–$(round(maximum(ao.ω̃), digits=4))"
                fᵗ_str = ao.fᵗ isa Number ? "$(round(ao.fᵗ, digits=4))" : "$(round(minimum(ao.fᵗ), digits=4))–$(round(maximum(ao.fᵗ), digits=4))"
                connector = (ia < length(aer_ops)) ? _TREE_MID : _TREE_END
                println(io, _TREE_PIPE, connector, "aerosol #$ia: ω̃=$ω_str, fᵗ=$fᵗ_str")
            end
        end
    else
        println(io, _TREE_PIPE, _TREE_END, "τ_rayl: $nBands band(s)")
    end

    # surfaces
    surf_strs = [sprint(show, s) for s in m.surfaces]
    println(io, _TREE_MID, "surfaces: ", join(surf_strs, ", "))

    # sources (v0.6 source-term refactor)
    print(io, _TREE_END, "sources: ", sprint(show, m.sources))
end

# Compact show (single-line, e.g. inside arrays)
function Base.show(io::IO, m::RTModel)
    FT = float_type(m)
    nBands = length(m.atmosphere.spec_bands)
    nSources = m.sources isa SourceSet ? length(m.sources) : (m.sources isa NoSource ? 0 : 1)
    print(io, "RTModel{$(typeof(m.architecture)), $FT}($nBands band(s), $nSources source(s))")
end

# ── vSmartMOM_Parameters ──────────────────────────────────────────────────

function Base.summary(io::IO, x::vSmartMOM_Parameters{FT}) where FT
    nBands = length(x.spec_bands)
    nLevels = length(x.p)
    print(io, "vSmartMOM_Parameters{$FT}($nBands band(s), $nLevels levels)")
end

function Base.show(io::IO, ::MIME"text/plain", x::vSmartMOM_Parameters{FT}) where FT
    nBands = length(x.spec_bands)

    println(io, "vSmartMOM_Parameters{$FT}")

    # ── radiative_transfer ──
    println(io, _TREE_MID, "radiative_transfer")
    # spectral bands
    band_strs = [begin
        b = x.spec_bands[i]
        "[$(round(b[1], digits=1))–$(round(b[end], digits=1))]"
    end for i in 1:min(nBands, 4)]
    band_summary = join(band_strs, " ")
    if nBands > 4
        band_summary *= " …"
    end
    println(io, _TREE_PIPE, _TREE_MID, "spectral_bands: $nBands band(s) $band_summary cm⁻¹")
    # surfaces
    surf_strs = [sprint(show, s) for s in x.brdf]
    println(io, _TREE_PIPE, _TREE_MID, "surfaces: ", join(surf_strs, ", "))
    # quadrature & polarization
    quad_name = typeof(x.quadrature_type).name.name
    pol_name = typeof(x.polarization_type).name.name
    println(io, _TREE_PIPE, _TREE_MID, "quadrature: $quad_name")
    println(io, _TREE_PIPE, _TREE_MID, "polarization: $pol_name")
    println(io, _TREE_PIPE, _TREE_END, "max_m=$(x.max_m), l_trunc=$(x.l_trunc), Δ_angle=$(x.Δ_angle)°, depol=$(x.depol)")

    # ── geometry ──
    vza_str = length(x.vza) <= 4 ? string(x.vza) : "[$(length(x.vza)) angles]"
    vaz_str = length(x.vaz) <= 4 ? string(x.vaz) : "[$(length(x.vaz)) angles]"
    println(io, _TREE_MID, "geometry: SZA=$(x.sza)°, VZA=$(vza_str)°, VAZ=$(vaz_str)°")

    # ── atmosphere ──
    nLevels = length(x.p)
    red_str = x.profile_reduction_n == -1 ? "" : ", reduced to $(x.profile_reduction_n) layers"
    println(io, _TREE_MID, "atmosphere: $nLevels levels$red_str")

    # ── absorption ──
    if !isnothing(x.absorption_params)
        ap = x.absorption_params
        mol_strs = unique(vcat(ap.fixed_molecules..., ap.variable_molecules...))
        h2o_active = !isempty(ap.h2o_lut) && any(!isnothing, ap.h2o_lut)
        all_mols = h2o_active ? vcat(["H2O (q-driven)"], mol_strs) : mol_strs
        bf_name = typeof(ap.broadening_function).name.name
        println(io, _TREE_MID, "absorption: $(length(all_mols)) molecule(s) [$(join(all_mols, ", "))], $bf_name, wing_cutoff=$(ap.wing_cutoff) cm⁻¹")
    else
        println(io, _TREE_MID, "absorption: none")
    end

    # ── scattering ──
    if !isnothing(x.scattering_params)
        sp = x.scattering_params
        nAer = length(sp.rt_aerosols)
        decomp_name = typeof(sp.decomp_type).name.name
        print(io, _TREE_END, "scattering: $nAer aerosol(s), λ_ref=$(sp.λ_ref) μm, $decomp_name")
        if nAer > 0
            println(io)
            for (i, rta) in enumerate(sp.rt_aerosols)
                connector = i < nAer ? _TREE_MID : _TREE_END
                dist = rta.aerosol.size_distribution
                dist_str = sprint(show, dist)
                print(io, _TREE_SPACE, connector,
                      "aerosol #$i: τ_ref=$(rta.τ_ref), $dist_str, nᵣ=$(rta.aerosol.nᵣ), nᵢ=$(rta.aerosol.nᵢ)")
                if i < nAer
                    println(io)
                end
            end
        end
    else
        print(io, _TREE_END, "scattering: none")
    end
end
