# Full-grid GCHP aerosol optical-depth diagnostic.
#
# This is intentionally outside CoreRT: it validates the microphysics → Mie
# extinction boundary without changing the production RT aerosol path.

const _AOD_M_AIR_KG_PER_MOL = 28.9644e-3

"""
    write_gchp_aod_diagnostic(out_path, gchp_path; wavelengths_um, ...)

Compute aerosol optical depth for every requested GCHP grid cell and write a
NetCDF file with dimensions `(x, y, face, wavelength)`.

The default wavelengths are `[0.760, 1.600, 2.200]` μm. The GCHP file is
opened once and aerosol ingest defaults to `:auto`.
"""
function write_gchp_aod_diagnostic(out_path::AbstractString,
                                   gchp_path::AbstractString;
                                   wavelengths_um = [0.760, 1.600, 2.200],
                                   aerosol_scheme = :auto,
                                   faces = nothing,
                                   xrange = nothing,
                                   yrange = nothing,
                                   integration = nothing,
                                   mie_lut = nothing,
                                   ri_round_digits = nothing,
                                   FT::DataType = Float64)
    λs = FT.(wavelengths_um)
    mie_cache = Dict{Any, FT}()

    GCHPFile(gchp_path; aerosol_scheme=aerosol_scheme, FT=FT) do f
        face_iter = faces === nothing ? (1:f.nf) : faces
        x_iter    = xrange === nothing ? (1:f.nx) : xrange
        y_iter    = yrange === nothing ? (1:f.ny) : yrange

        NCDataset(String(out_path), "c") do ds
            defDim(ds, "x", f.nx)
            defDim(ds, "y", f.ny)
            defDim(ds, "face", f.nf)
            defDim(ds, "wavelength", length(λs))

            λ_var = defVar(ds, "wavelength_um", Float64, ("wavelength",))
            λ_var[:] = Float64.(λs)
            λ_var.attrib["units"] = "micrometer"
            λ_var.attrib["long_name"] = "wavelength"

            lat = defVar(ds, "lat", Float64, ("x", "y", "face"))
            lon = defVar(ds, "lon", Float64, ("x", "y", "face"))
            lat[:, :, :] = f.lats
            lon[:, :, :] = f.lons
            lat.attrib["units"] = "degree_north"
            lon.attrib["units"] = "degree_east"

            aod = defVar(ds, "aod", Float64, ("x", "y", "face", "wavelength"))
            aod[:, :, :, :] .= NaN
            aod.attrib["long_name"] = "aerosol optical depth from sectional aerosol microphysics"
            aod.attrib["units"] = "1"

            for idf in face_iter, idy in y_iter, idx in x_iter
                scene = scene_at(f, idx, idy, idf; FT=FT)
                aod_kwargs = integration === nothing ?
                    (; mie_cache=mie_cache, mie_lut=mie_lut,
                     ri_round_digits=ri_round_digits) :
                    (; integration=integration, mie_cache=mie_cache,
                     mie_lut=mie_lut, ri_round_digits=ri_round_digits)
                aod[idx, idy, idf, :] = Float64.(compute_scene_aod(
                    scene, λs; aod_kwargs...))
            end

            ds.attrib["source_file"] = String(gchp_path)
            ds.attrib["aerosol_scheme"] = string(aerosol_scheme)
            ds.attrib["time"] = string(f.time)
            ds.attrib["aod_algorithm"] = "Mie extinction at sectional bin centers with volume-weighted effective RI"
            ds.attrib["ri_round_digits"] = ri_round_digits === nothing ? "none" : string(ri_round_digits)
            ds.attrib["mie_lut"] = mie_lut === nothing ? "none" : string(typeof(mie_lut))
            ds.attrib["integration"] = integration === nothing ? "scene default" : string(integration)
        end
    end

    return String(out_path)
end

function _range_chunks(r, chunk_size::Integer)
    vals = collect(r)
    isempty(vals) && return UnitRange{Int}[]
    chunks = UnitRange{Int}[]
    i = firstindex(vals)
    while i <= lastindex(vals)
        j = min(i + Int(chunk_size) - 1, lastindex(vals))
        push!(chunks, vals[i]:vals[j])
        i = j + 1
    end
    return chunks
end

function _bulk_var3(ds, name::AbstractString, xr, yr, idf, FT)
    return FT.(Array(ds[name].var[xr, yr, idf, :, 1]))
end

function _bulk_var2(ds, name::AbstractString, xr, yr, idf, FT)
    return FT.(Array(ds[name].var[xr, yr, idf, 1]))
end

function _chunk_sectional_data(scheme::TOMASScheme{FT}, ds, xr, yr, idf;
                               integration,
                               mixing_rule = VolumeWeightedMixing()) where {FT}
    n_bins = scheme.n_bins
    nlev = length(ds["lev"])
    species = scheme.species
    nspc = length(species)

    Met_AD = _bulk_var3(ds, "Met_AD", xr, yr, idf, FT)
    Met_AIRVOL = _bulk_var3(ds, "Met_AIRVOL", xr, yr, idf, FT)
    air = Met_AD ./ FT(_AOD_M_AIR_KG_PER_MOL)
    nx = length(xr)
    ny = length(yr)

    number = zeros(FT, n_bins, nlev, nx, ny)
    for ibin in 1:n_bins
        var_name = @sprintf("SpeciesConcVV_NK%02d", ibin)
        haskey(ds, var_name) || continue
        nk = _bulk_var3(ds, var_name, xr, yr, idf, FT)
        @inbounds for iy in 1:ny, ix in 1:nx, ilev in 1:nlev
            itoa = nlev - ilev + 1
            number[ibin, itoa, ix, iy] =
                (nk[ix, iy, ilev] / FT(1000)) * air[ix, iy, ilev] /
                (Met_AIRVOL[ix, iy, ilev] * FT(1e6))
        end
    end

    mass = zeros(FT, n_bins, nlev, nspc, nx, ny)
    for (ispc, spc) in enumerate(species)
        MW_ug_per_mol = scheme.molar_masses[spc] * FT(1e9)
        for ibin in 1:n_bins
            var_name = string("SpeciesConcVV_", spc, @sprintf("%02d", ibin))
            haskey(ds, var_name) || continue
            vv = _bulk_var3(ds, var_name, xr, yr, idf, FT)
            @inbounds for iy in 1:ny, ix in 1:nx, ilev in 1:nlev
                itoa = nlev - ilev + 1
                mass[ibin, itoa, ispc, ix, iy] =
                    vv[ix, iy, ilev] * air[ix, iy, ilev] * MW_ug_per_mol /
                    Met_AIRVOL[ix, iy, ilev]
            end
        end
    end

    return (; number, mass)
end

function _bulk_profile_chunk(ds, xr, yr, idf, FT)
    dp = _bulk_var3(ds, "Met_DELP", xr, yr, idf, FT)
    sp = _bulk_var2(ds, "Met_PS2WET", xr, yr, idf, FT)
    T_boa = _bulk_var3(ds, "Met_T", xr, yr, idf, FT)
    q_raw = _bulk_var3(ds, "Met_SPHU", xr, yr, idf, FT)
    q_units = ds["Met_SPHU"].attrib["units"]
    q_scale = q_units == "g kg-1" ? FT(1e-3) :
              q_units == "kg kg-1" ? one(FT) :
              throw(ArgumentError("Unsupported Met_SPHU units: '$q_units'"))

    return (; dp, sp, T_boa, q_raw, q_scale)
end

function _column_p_half(dp_chunk, sp_chunk, ix::Integer, iy::Integer, FT)
    dp_col = view(dp_chunk, ix, iy, :)
    sp = sp_chunk[ix, iy]
    return FT.(reverse([sp; sp .+ cumsum(-collect(dp_col))]))
end

function _aod_size_to_um(value::FT, units::AbstractString) where {FT}
    if units in ("μm", "um", "micrometer", "micrometers")
        return value
    elseif units in ("nm", "nanometer", "nanometers")
        return value / FT(1000)
    elseif units in ("m", "meter", "meters")
        return value * FT(1e6)
    end
    throw(ArgumentError("Cannot convert aerosol size-grid units '$units' to μm"))
end

function _aod_bin_log_bounds(grid, ibin::Integer, FT)
    grid.edges === nothing &&
        throw(ArgumentError("bulk AOD bin integration requires explicit size-grid edges"))
    lo = _aod_size_to_um(FT(grid.edges[ibin]), grid.units)
    hi = _aod_size_to_um(FT(grid.edges[ibin + 1]), grid.units)
    lo > zero(FT) && hi > lo ||
        throw(ArgumentError("invalid aerosol bin edges for bin $ibin"))
    return log10(lo), log10(hi)
end

function _aod_radius_from_log_size(grid, log_size::FT) where {FT}
    value_um = FT(10)^log_size
    if grid.coordinate === :diameter
        return value_um / FT(2)
    elseif grid.coordinate === :radius
        return value_um
    end
    throw(ArgumentError("AOD diagnostic needs radius/diameter size grid; got coordinate :$(grid.coordinate)"))
end

function _aod_radius_center(grid, ibin::Integer, FT)
    value_um = _aod_size_to_um(FT(grid.centers[ibin]), grid.units)
    if grid.coordinate === :diameter
        return value_um / FT(2)
    elseif grid.coordinate === :radius
        return value_um
    end
    throw(ArgumentError("AOD diagnostic needs radius/diameter size grid; got coordinate :$(grid.coordinate)"))
end

function _aod_mean_dndlog(number_col, grid, ibin::Integer, ilev::Integer, FT)
    xlo, xhi = _aod_bin_log_bounds(grid, ibin, FT)
    width = xhi - xlo
    width > zero(FT) || return zero(FT)
    return FT(number_col[ibin, ilev]) * FT(1e6) / width
end

function _aod_limited_linear_slope(number_col, grid, ibin::Integer, ilev::Integer,
                                   mean_i::FT, xlo::FT, xmid::FT,
                                   xhi::FT) where {FT}
    n_bins = size(number_col, 1)
    if ibin == 1
        mean_i <= zero(FT) && return zero(FT)
        return mean_i / (xmid - xlo)
    elseif ibin == n_bins
        mean_i <= zero(FT) && return zero(FT)
        return -mean_i / (xhi - xmid)
    end

    xprev = sum(_aod_bin_log_bounds(grid, ibin - 1, FT)) / FT(2)
    xnext = sum(_aod_bin_log_bounds(grid, ibin + 1, FT)) / FT(2)
    mean_prev = _aod_mean_dndlog(number_col, grid, ibin - 1, ilev, FT)
    mean_next = _aod_mean_dndlog(number_col, grid, ibin + 1, ilev, FT)
    s_left = (mean_i - mean_prev) / (xmid - xprev)
    s_right = (mean_next - mean_i) / (xnext - xmid)

    slope = s_left * s_right <= zero(FT) ?
        zero(FT) :
        sign(s_left) * min(abs(s_left), abs(s_right))

    if mean_i > zero(FT)
        max_abs_slope = mean_i / max(abs(xlo - xmid), abs(xhi - xmid))
        slope = clamp(slope, -max_abs_slope, max_abs_slope)
    else
        slope = zero(FT)
    end
    return slope
end

function _bulk_integrated_bin_extinction(number_col, grid, ibin::Integer,
                                         ilev::Integer, λ::FT,
                                         n_eff::Complex{FT},
                                         ::DirectBinSum, mie_cache,
                                         mie_lut, ri_round_digits) where {FT}
    N_m3 = FT(number_col[ibin, ilev]) * FT(1e6)
    r_um = _aod_radius_center(grid, ibin, FT)
    c_ext = Aerosols._cached_mie_extinction(r_um, λ, n_eff, mie_cache,
                                            mie_lut, ri_round_digits)
    return N_m3 * c_ext
end

function _bulk_integrated_bin_extinction(number_col, grid, ibin::Integer,
                                         ilev::Integer, λ::FT,
                                         n_eff::Complex{FT},
                                         integration::ConstantIntegrationPerBin,
                                         mie_cache, mie_lut,
                                         ri_round_digits) where {FT}
    xlo, xhi = _aod_bin_log_bounds(grid, ibin, FT)
    xmid = (xlo + xhi) / FT(2)
    half_width = (xhi - xlo) / FT(2)
    mean_dndlog = _aod_mean_dndlog(number_col, grid, ibin, ilev, FT)
    nodes, weights = Aerosols._gausslegendre_nodes_weights(integration.nquad)

    β = zero(FT)
    @inbounds for i in eachindex(nodes)
        x = xmid + half_width * FT(nodes[i])
        r_um = _aod_radius_from_log_size(grid, x)
        c_ext = Aerosols._cached_mie_extinction(r_um, λ, n_eff, mie_cache,
                                                mie_lut, ri_round_digits)
        β += FT(weights[i]) * mean_dndlog * c_ext
    end
    return half_width * β
end

function _bulk_integrated_bin_extinction(number_col, grid, ibin::Integer,
                                         ilev::Integer, λ::FT,
                                         n_eff::Complex{FT},
                                         integration::LinearIntegrationPerBin,
                                         mie_cache, mie_lut,
                                         ri_round_digits) where {FT}
    xlo, xhi = _aod_bin_log_bounds(grid, ibin, FT)
    xmid = (xlo + xhi) / FT(2)
    half_width = (xhi - xlo) / FT(2)
    mean_dndlog = _aod_mean_dndlog(number_col, grid, ibin, ilev, FT)
    slope = _aod_limited_linear_slope(number_col, grid, ibin, ilev,
                                      mean_dndlog, xlo, xmid, xhi)
    nodes, weights = Aerosols._gausslegendre_nodes_weights(integration.nquad)

    β = zero(FT)
    @inbounds for i in eachindex(nodes)
        x = xmid + half_width * FT(nodes[i])
        dndlog = max(zero(FT), mean_dndlog + slope * (x - xmid))
        r_um = _aod_radius_from_log_size(grid, x)
        c_ext = Aerosols._cached_mie_extinction(r_um, λ, n_eff, mie_cache,
                                                mie_lut, ri_round_digits)
        β += FT(weights[i]) * dndlog * c_ext
    end
    return half_width * β
end

function _bulk_column_aod(number_col, mass_col, scheme::TOMASScheme{FT},
                          p_half, T, q, λs;
                          integration,
                          mie_cache,
                          mie_lut,
                          ri_round_digits,
                          ri_speciesλ,
                          density_by_species) where {FT}
    nlev = size(number_col, 2)
    out = zeros(FT, length(λs))
    dz = [Aerosols._layer_thickness_m(FT(p_half[ilev]), FT(p_half[ilev + 1]),
                                      FT(T[ilev]), FT(q[ilev]))
          for ilev in 1:nlev]

    nspc = length(scheme.species)
    vol = zeros(FT, nspc)
    @inbounds for (iλ, λ) in enumerate(λs)
        τ = zero(FT)
        for ilev in 1:nlev
            dz_layer = dz[ilev]
            dz_layer > zero(FT) || continue
            for ibin in 1:scheme.n_bins
                number_col[ibin, ilev] > zero(FT) || continue
                total_vol = zero(FT)
                for ispc in 1:nspc
                    v = FT(mass_col[ibin, ilev, ispc]) / density_by_species[ispc]
                    vol[ispc] = v
                    total_vol += v
                end
                total_vol > zero(FT) || continue

                n_eff = Complex{FT}(zero(FT), zero(FT))
                for ispc in 1:nspc
                    n_eff += (vol[ispc] / total_vol) * ri_speciesλ[ispc, iλ]
                end
                τ += _bulk_integrated_bin_extinction(
                    number_col, scheme.size_grid, ibin, ilev, λ, n_eff,
                    integration, mie_cache, mie_lut, ri_round_digits) * dz_layer
            end
        end
        out[iλ] = max(zero(FT), τ)
    end
    return out
end

function _ri_species_wavelengths(scheme::TOMASScheme{FT}, λs) where {FT}
    db = Aerosols._default_refractive_index_database(FT)
    out = Matrix{Complex{FT}}(undef, length(scheme.species), length(λs))
    for (ispc, spc) in enumerate(scheme.species), (iλ, λ) in enumerate(λs)
        key = refractive_index_key(scheme, spc)
        out[ispc, iλ] = get_refractive_index(db, key, λ)
    end
    return out
end

"""
    write_gchp_aod_diagnostic_bulk(out_path, gchp_path; ...)

Chunked version of [`write_gchp_aod_diagnostic`](@ref). It reads each requested
GCHP face/chunk variable as an array, then computes all columns in memory. This
avoids the per-cell NetCDF read pattern in `scene_at` and is intended for fast
global AOD diagnostics.

Currently this bulk path supports TOMAS-family sectional schemes.
"""
function write_gchp_aod_diagnostic_bulk(out_path::AbstractString,
                                        gchp_path::AbstractString;
                                        wavelengths_um = [0.760, 1.600, 2.200],
                                        aerosol_scheme = :auto,
                                        faces = nothing,
                                        xrange = nothing,
                                        yrange = nothing,
                                        xchunk::Integer = 15,
                                        ychunk::Integer = 90,
                                        integration = LinearIntegrationPerBin(5),
                                        mie_lut = nothing,
                                        ri_round_digits = nothing,
                                        threaded::Bool = true,
                                        FT::DataType = Float64)
    λs = FT.(wavelengths_um)
    mie_cache = Dict{Any, FT}()

    GCHPFile(gchp_path; aerosol_scheme=aerosol_scheme, FT=FT) do f
        scheme = f.aer_scheme
        scheme isa TOMASScheme ||
            throw(ArgumentError("bulk GCHP AOD currently requires a TOMASScheme; got $(typeof(scheme))"))

        face_iter = faces === nothing ? (1:f.nf) : faces
        x_iter = xrange === nothing ? (1:f.nx) : xrange
        y_iter = yrange === nothing ? (1:f.ny) : yrange

        x_chunks = _range_chunks(x_iter, xchunk)
        y_chunks = _range_chunks(y_iter, ychunk)
        ri_speciesλ = _ri_species_wavelengths(scheme, λs)
        density_by_species = densities(scheme, scheme.species)
        if integration isa Union{ConstantIntegrationPerBin, LinearIntegrationPerBin}
            Aerosols._gausslegendre_nodes_weights(integration.nquad)
        end

        NCDataset(String(out_path), "c") do ods
            defDim(ods, "x", f.nx)
            defDim(ods, "y", f.ny)
            defDim(ods, "face", f.nf)
            defDim(ods, "wavelength", length(λs))

            λ_var = defVar(ods, "wavelength_um", Float64, ("wavelength",))
            λ_var[:] = Float64.(λs)
            λ_var.attrib["units"] = "micrometer"
            λ_var.attrib["long_name"] = "wavelength"

            lat = defVar(ods, "lat", Float64, ("x", "y", "face"))
            lon = defVar(ods, "lon", Float64, ("x", "y", "face"))
            lat[:, :, :] = f.lats
            lon[:, :, :] = f.lons
            lat.attrib["units"] = "degree_north"
            lon.attrib["units"] = "degree_east"

            aod = defVar(ods, "aod", Float64, ("x", "y", "face", "wavelength"))
            aod[:, :, :, :] .= NaN
            aod.attrib["long_name"] = "aerosol optical depth from sectional aerosol microphysics"
            aod.attrib["units"] = "1"

            for idf in face_iter, yr in y_chunks, xr in x_chunks
                profile = _bulk_profile_chunk(f.ds, xr, yr, idf, FT)
                aero_chunk = _chunk_sectional_data(scheme, f.ds, xr, yr, idf;
                                                   integration=integration)
                nx = length(xr)
                ny = length(yr)
                aod_chunk = Array{FT}(undef, nx, ny, length(λs))

                function compute_column!(linear_index)
                    ix = (linear_index - 1) % nx + 1
                    iy = (linear_index - 1) ÷ nx + 1
                    p_half = _column_p_half(profile.dp, profile.sp, ix, iy, FT)
                    T = FT.(reverse(view(profile.T_boa, ix, iy, :)))
                    q = FT.(reverse(view(profile.q_raw, ix, iy, :))) .* profile.q_scale
                    local_cache = mie_lut === nothing ? Dict{Any, FT}() : mie_cache
                    aod_chunk[ix, iy, :] .= _bulk_column_aod(
                        view(aero_chunk.number, :, :, ix, iy),
                        view(aero_chunk.mass, :, :, :, ix, iy),
                        scheme, p_half, T, q, λs;
                        integration=integration, mie_cache=local_cache,
                        mie_lut=mie_lut, ri_round_digits=ri_round_digits,
                        ri_speciesλ=ri_speciesλ,
                        density_by_species=density_by_species)
                    return nothing
                end

                if threaded && Threads.nthreads() > 1
                    Threads.@threads for linear_index in 1:(nx * ny)
                        compute_column!(linear_index)
                    end
                else
                    for linear_index in 1:(nx * ny)
                        compute_column!(linear_index)
                    end
                end

                aod[xr, yr, idf, :] = Float64.(aod_chunk)
            end

            ods.attrib["source_file"] = String(gchp_path)
            ods.attrib["aerosol_scheme"] = string(aerosol_scheme)
            ods.attrib["time"] = string(f.time)
            ods.attrib["aod_algorithm"] = "bulk NetCDF chunk TOMAS Mie extinction with volume-weighted effective RI"
            ods.attrib["ri_round_digits"] = ri_round_digits === nothing ? "none" : string(ri_round_digits)
            ods.attrib["mie_lut"] = mie_lut === nothing ? "none" : string(typeof(mie_lut))
            ods.attrib["integration"] = string(integration)
            ods.attrib["xchunk"] = string(xchunk)
            ods.attrib["ychunk"] = string(ychunk)
            ods.attrib["threaded"] = string(threaded)
        end
    end

    return String(out_path)
end
