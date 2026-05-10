# =============================================================================
# mtckd.jl — MT_CKD water-vapor continuum (self + foreign).
#
# Reads the AER-distributed reference NetCDF (`absco-ref_wv-mt-ckd.nc`,
# MT_CKD v4.2/4.3) and adds the H₂O continuum optical depth to τ_abs at
# model build time. Implementation follows the standard LBLRTM convention:
#
#   σ_self(ν, T, p_h2o) = C_self_ref(ν) · radterm(ν, T) ·
#                          (p_h2o / p_ref) · (T_ref / T)^texp(ν)
#   σ_for(ν, T, p_dry)  = C_for_ref(ν)  · radterm(ν, T) · (p_dry / p_ref)
#   τ_layer(ν)          = (σ_self + σ_for) · n_h2o · Δz
#
# where  radterm(ν, T) = ν · tanh(c₂ ν / 2T)  with c₂ = 1.4388 cm·K
# converts the AER "continuum coefficient" (cm²/molec/cm⁻¹) into a true
# cross-section (cm²/molec) — see the AER absco-ref Notes attribute and
# Mlawer et al. (2012), doi:10.1098/rsta.2011.0295.
#
# The AER table covers ν ∈ [-20, 20000] cm⁻¹ (≈ ∞–500 nm). Outside that
# range τ_continuum = 0 — matches the published model's domain.
# =============================================================================

const _MTCKD_C2 = 1.4388  # cm·K, h·c/k_B (radiation constant 2)

"""
    MTCKDTable

Water-vapor continuum reference table from MT_CKD's `absco-ref_wv-mt-ckd.nc`.
Stored at the file's native ν grid (uniform 10 cm⁻¹); per-band tables are
built on demand by `build_mtckd_band_table`.
"""
struct MTCKDTable
    ν::Vector{Float64}            # cm⁻¹, ascending
    C_self::Vector{Float64}       # cm²/molec/cm⁻¹ at T_ref
    C_for::Vector{Float64}        # cm²/molec/cm⁻¹ at T_ref
    self_texp::Vector{Float64}    # dimensionless
    p_ref::Float64                # mbar (= hPa)
    T_ref::Float64                # K
end

"""
    MTCKDBand

Per-band cache of MT_CKD coefficients pre-interpolated onto a model ν grid.
Filling values outside the file's ν range with zeros.
"""
struct MTCKDBand
    C_self::Vector{Float64}    # length(ν_grid)
    C_for::Vector{Float64}
    texp::Vector{Float64}
    p_ref::Float64
    T_ref::Float64
end

"""
    load_mtckd(path) -> MTCKDTable

Read the AER-distributed `absco-ref_wv-mt-ckd.nc` file.
"""
function load_mtckd(path::AbstractString)
    ν       = Float64.(NetCDF.ncread(path, "wavenumbers"))
    C_self  = Float64.(NetCDF.ncread(path, "self_absco_ref"))
    C_for   = Float64.(NetCDF.ncread(path, "for_absco_ref"))
    texp    = Float64.(NetCDF.ncread(path, "self_texp"))
    p_ref   = Float64(NetCDF.ncread(path, "ref_press")[])
    T_ref   = Float64(NetCDF.ncread(path, "ref_temp")[])
    return MTCKDTable(ν, C_self, C_for, texp, p_ref, T_ref)
end

"""
    build_mtckd_band(table, ν_grid) -> MTCKDBand

Interpolate the table's coefficients onto `ν_grid`. Outside the table's
ν range (typically [-20, 20000] cm⁻¹), all coefficients are set to 0 so
τ_continuum vanishes there (UV/Vis above 500 nm has no MT_CKD contribution).
"""
function build_mtckd_band(table::MTCKDTable, ν_grid::AbstractVector)
    nν = length(ν_grid)
    Cs = zeros(Float64, nν)
    Cf = zeros(Float64, nν)
    te = zeros(Float64, nν)
    νlo, νhi = table.ν[1], table.ν[end]
    @inbounds for (k, νq_) in enumerate(ν_grid)
        νq = Float64(νq_)
        (νq < νlo || νq > νhi) && continue
        j = searchsortedfirst(table.ν, νq)
        if j ≤ 1
            Cs[k] = table.C_self[1];     Cf[k] = table.C_for[1];    te[k] = table.self_texp[1]
        elseif j > length(table.ν)
            Cs[k] = table.C_self[end];   Cf[k] = table.C_for[end];  te[k] = table.self_texp[end]
        else
            ν1, ν2 = table.ν[j-1], table.ν[j]
            w = (νq - ν1) / (ν2 - ν1)
            Cs[k] = (1 - w) * table.C_self[j-1]    + w * table.C_self[j]
            Cf[k] = (1 - w) * table.C_for[j-1]     + w * table.C_for[j]
            te[k] = (1 - w) * table.self_texp[j-1] + w * table.self_texp[j]
        end
    end
    return MTCKDBand(Cs, Cf, te, table.p_ref, table.T_ref)
end

"""
    compute_τ_h2o_continuum!(τ_abs, band, ν_grid, profile, vmr_h2o)

Add the MT_CKD H₂O self+foreign continuum optical depth to `τ_abs[ν, layer]`.

`vmr_h2o` is per-layer or scalar. `profile.p_full` is in hPa, `profile.T`
in K, `profile.Δz` in m. All intermediate arithmetic is Float64; the
result is converted to `eltype(τ_abs)` only at accumulation.
"""
function compute_τ_h2o_continuum!(τ_abs::AbstractMatrix,
                                   band::MTCKDBand,
                                   ν_grid::AbstractVector,
                                   profile,
                                   vmr_h2o)
    nλ, nlay = size(τ_abs)
    nλ == length(band.C_self) ||
        error("compute_τ_h2o_continuum!: ν-grid length mismatch — band $(length(band.C_self)), τ_abs $(nλ)")
    FT_τ = eltype(τ_abs)
    p_ref = band.p_ref
    T_ref = band.T_ref
    @inbounds for iz in 1:nlay
        T   = Float64(profile.T[iz])
        P   = Float64(profile.p_full[iz])
        v_h = Float64(vmr_h2o isa AbstractVector ? vmr_h2o[iz] : vmr_h2o)
        # Number densities (molec/cm³)
        n_air = P * 1e2 / (1.380649e-23 * T) * 1e-6
        n_h2o = v_h * n_air
        # Partial pressures (hPa)
        p_h2o = v_h * P
        p_dry = P - p_h2o
        Δz_cm = Float64(profile.Δz[iz]) * 100.0
        # Layer-constant pressure scaling
        self_pscale = p_h2o / p_ref
        for_pscale  = p_dry / p_ref
        # Loop over ν
        for k in 1:nλ
            ν       = Float64(ν_grid[k])
            radterm = ν * tanh(_MTCKD_C2 * ν / (2.0 * T))
            σ_self  = band.C_self[k] * radterm * self_pscale *
                      (T_ref / T)^band.texp[k]
            σ_for   = band.C_for[k]  * radterm * for_pscale
            τ_abs[k, iz] += FT_τ((σ_self + σ_for) * n_h2o * Δz_cm)
        end
    end
    return τ_abs
end

"""
    compute_τ_h2o_continuum!(τ_abs, table::MTCKDTable, ν_grid, profile, vmr_h2o)

Convenience overload that builds the per-band cache on the fly.
"""
function compute_τ_h2o_continuum!(τ_abs::AbstractMatrix,
                                   table::MTCKDTable,
                                   ν_grid::AbstractVector,
                                   profile,
                                   vmr_h2o)
    band = build_mtckd_band(table, ν_grid)
    return compute_τ_h2o_continuum!(τ_abs, band, ν_grid, profile, vmr_h2o)
end
