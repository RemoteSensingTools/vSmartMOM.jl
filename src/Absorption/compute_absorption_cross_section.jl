#=

This file contains functions that perform the core absorption cross section calculations.
`compute_absorption_cross_section` is the main entry point, also wrapped by
`absorption_cross_section` in autodiff_helper.jl.

`compute_absorption_cross_section` is implemented both from scratch (::HitranModel) and as
an interpolation (::InterpolationModel). For the HitranModel, a batched kernel processes
all spectral lines in a single GPU/CPU launch for efficiency.

=#

"""
    compute_absorption_cross_section(model::HitranModel, grid, pressure, temperature; wavelength_flag=false)

Compute absorption cross-section from HITRAN line-by-line data.

Sums Voigt (or Doppler/Lorentz) line shapes over all transitions in the HITRAN table
within the grid range. Uses a batched kernel to process all lines in a single launch
for efficient GPU/CPU execution.

# Arguments
- `model::HitranModel`: Model with HITRAN data, broadening, wing_cutoff, CEF
- `grid`: Wavenumber [cm⁻¹] or wavelength [nm] grid
- `pressure::Real`: Pressure (hPa)
- `temperature::Real`: Temperature (K)
- `wavelength_flag::Bool=false`: If true, grid is wavelength in nm

# Returns
- `Array`: Absorption cross-section (cm²/molecule) at each grid point
"""
function compute_absorption_cross_section(
                # Required
                model::HitranModel,          # Model to use in this cross section calculation
                                             # (Calculation from Hitran data vs. using Interpolator)
                # Wavelength [nm] or wavenumber [cm-1] grid (modify using wavelength_flag)
                grid::Union{AbstractRange{<:Real}, AbstractArray},  # Can be range OR array
                pressure::Real,              # actual pressure [hPa]
                temperature::Real;           # actual temperature [K]
                # Optionals
                wavelength_flag::Bool=false  # Use wavelength in nm (true) or wavenumber cm-1 units (false)
                )

    (; hitran, broadening, wing_cutoff, vmr, CEF, architecture) = model

    # Notify user of wavelength grid
    if (wavelength_flag)
        @info """
        Note: Cross-section reported to wavelength grid
        (will internally convert to wavenumber grid for calculations)
        """  maxlog = 5
    end

    # Convert T to float type (ex. if Int)
    temperature = AbstractFloat(temperature)

    FT = eltype(temperature)

    # Store results here to return
    result = array_type(architecture)(zeros(FT, length(grid)))

    # Calculate the minimum and maximum grid bounds, including the wing cutoff
    grid_max = maximum(grid) + wing_cutoff
    grid_min = minimum(grid) - wing_cutoff

    # Convert to wavenumber from [nm] space if necessary
    grid = wavelength_flag ? reverse(nm_per_m ./ grid) : grid
    grid_min, grid_max = wavelength_flag ? (nm_per_m /grid_max, nm_per_m/grid_min) : (grid_min, grid_max)

    # Interpolators from grid bounds to index values
    N_grid = length(grid)
    if N_grid > 1
        grid_idx_interp_low  = LinearInterpolation(grid, 1:1:N_grid, extrapolation_bc=1)
        grid_idx_interp_high = LinearInterpolation(grid, 1:1:N_grid, extrapolation_bc=N_grid)
    end

    # --- Pre-compute line parameters on CPU (Steps 1b + 3: batch + cache qoft) ---

    # Cache partition sum ratios per unique (mol, iso) pair
    qoft_cache = Dict{Tuple{Int,Int}, FT}()
    rate = zeros(FT, 1)

    # Count active lines first for pre-allocation
    n_active = 0
    for j in eachindex(hitran.Sᵢ)
        if grid_min < hitran.νᵢ[j] < grid_max
            n_active += 1
        end
    end

    # Pre-allocate parameter arrays
    ν_arr       = Vector{FT}(undef, n_active)
    γ_d_arr     = Vector{FT}(undef, n_active)
    y_arr       = Vector{FT}(undef, n_active)
    S_arr       = Vector{FT}(undef, n_active)
    γ_l_arr     = Vector{FT}(undef, n_active)
    istart_arr  = Vector{Int32}(undef, n_active)
    istop_arr   = Vector{Int32}(undef, n_active)

    # Fill parameter arrays
    k = 0
    for j in eachindex(hitran.Sᵢ)
        if !(grid_min < hitran.νᵢ[j] < grid_max)
            continue
        end
        k += 1

        # Apply pressure shift
        ν   = hitran.νᵢ[j] + pressure / p_ref * hitran.δ_air[j]

        # Compute Lorentzian HWHM
        γ_l = (hitran.γ_air[j] *
              (1 - vmr) * pressure / p_ref + hitran.γ_self[j] *
              vmr * pressure / p_ref) *
              (t_ref / temperature)^hitran.n_air[j]

        # Compute Doppler HWHM
        γ_d = ((cSqrt2Ln2 / cc_) * sqrt(cBolts_ / cMassMol) * sqrt(temperature) *
        hitran.νᵢ[j] / sqrt(mol_weight(hitran.mol[j], hitran.iso[j])))

        # Ratio of widths
        y = sqrt(cLn2) * γ_l / γ_d

        # Apply line intensity temperature corrections
        S = hitran.Sᵢ[j]
        if hitran.E″[j] != -1
            # Cache qoft result per unique (mol, iso) pair
            key = (hitran.mol[j], hitran.iso[j])
            if !haskey(qoft_cache, key)
                qoft!(hitran.mol[j], hitran.iso[j], temperature, t_ref, rate)
                qoft_cache[key] = rate[1]
            end
            S = S * qoft_cache[key] *
                    exp(c₂ * hitran.E″[j] * (1 / t_ref - 1 / temperature)) *
                    (-expm1(-c₂ * hitran.νᵢ[j] / temperature)) / (-expm1(-c₂ * hitran.νᵢ[j] / t_ref))
        end

        # Calculate index range that this line impacts
        if N_grid > 1
            istart = Int32(round(grid_idx_interp_low(ν - wing_cutoff)))
            istop  = Int32(round(grid_idx_interp_high(ν + wing_cutoff)))
        else
            istart = Int32(1)
            istop  = Int32(1)
        end

        ν_arr[k]      = ν
        γ_d_arr[k]    = γ_d
        y_arr[k]      = y
        S_arr[k]      = S
        γ_l_arr[k]    = γ_l
        istart_arr[k] = istart
        istop_arr[k]  = istop
    end

    # --- Launch batched kernel ---
    if n_active > 0 && N_grid > 0
        device = devi(architecture)

        # Upload grid and parameter arrays to device
        grid_d    = array_type(architecture)(grid)
        ν_d       = array_type(architecture)(ν_arr)
        γ_d_d     = array_type(architecture)(γ_d_arr)
        y_d       = array_type(architecture)(y_arr)
        S_d       = array_type(architecture)(S_arr)
        γ_l_d     = array_type(architecture)(γ_l_arr)
        istart_d  = array_type(architecture)(istart_arr)
        istop_d   = array_type(architecture)(istop_arr)

        kernel! = line_shape_batch!(device)
        kernel!(result, grid_d, ν_d, γ_d_d, γ_l_d, y_d, S_d,
                istart_d, istop_d, Int32(n_active), broadening, CEF,
                ndrange=N_grid)
        synchronize_if_gpu()
    end

    # Return the resulting lineshape
    return (wavelength_flag ? reverse(result) : result)
end

"""
    compute_absorption_cross_section(model::InterpolationModel, grid, pressure, temperature; wavelength_flag=false)

Compute absorption cross-section by interpolating a pre-computed lookup table.

Uses BSpline interpolation over (ν, p, T). Faster than line-by-line for repeated calls.

# Arguments
- `model::InterpolationModel`: Pre-computed interpolator with ν_grid, p_grid, t_grid
- `grid`: Wavenumber [cm⁻¹] or wavelength [nm] grid
- `pressure::Real`: Pressure (hPa)
- `temperature::Real`: Temperature (K)
- `wavelength_flag::Bool=false`: If true, grid is wavelength in nm

# Returns
- `Array`: Interpolated absorption cross-section (cm²/molecule) at each grid point
"""
function compute_absorption_cross_section(
    # Required
    model::InterpolationModel,      # Model to use in this cross section calculation
                                    # (Calculation from Interpolator vs. Hitran Data)
    # Wavelength [nm] or wavenumber [cm-1] grid
    grid::Union{AbstractRange{<:Real}, AbstractArray},  # Can be range OR array
    pressure::Real,                 # actual pressure [hPa]
    temperature::Real;              # actual temperature [K]
    # Optionals
    wavelength_flag::Bool=false,    # Use wavelength in nm (true) or wavenumber cm-1 units (false)
    )

    # Convert to wavenumber from [nm] space if necessary
    grid = wavelength_flag ? reverse(nm_per_m ./ collect(grid)) : collect(grid)

    # Scale the interpolation to match the model grids
    sitp = scale(model.itp, model.ν_grid, model.p_grid, model.t_grid)

    # Perform the interpolation and return the resulting grid
    return sitp(grid, pressure, temperature)
end

#=

Batched line-shape kernels. Each thread handles one grid point and loops over all
active spectral lines, accumulating contributions. This replaces the per-line kernel
launch pattern with a single kernel launch.

=#

"""
    line_shape_batch!(A, grid, ν_arr, γ_d_arr, γ_l_arr, y_arr, S_arr,
                      istart_arr, istop_arr, N_lines, ::Doppler, CEF)

KernelAbstractions batched Doppler line-shape kernel. Each workitem owns one
spectral grid index, loops over all active HITRAN lines, applies each line's
wing-cutoff index window, and accumulates the Gaussian contribution into
`A[I]`.
"""
@kernel function line_shape_batch!(A, @Const(grid), @Const(ν_arr), @Const(γ_d_arr),
                                   @Const(γ_l_arr), @Const(y_arr), @Const(S_arr),
                                   @Const(istart_arr), @Const(istop_arr),
                                   N_lines, ::Doppler, CEF)
    I = @index(Global, Linear)
    FT = eltype(A)
    acc = zero(FT)
    ν_i = FT(grid[I])
    @inbounds for j in 1:N_lines
        if istart_arr[j] <= I <= istop_arr[j]
            acc += FT(S_arr[j]) * FT(cSqrtLn2divSqrtPi) *
                   exp(-FT(cLn2) * ((ν_i - FT(ν_arr[j])) / FT(γ_d_arr[j]))^2) / FT(γ_d_arr[j])
        end
    end
    @inbounds A[I] += acc
end

"""
    line_shape_batch!(A, grid, ν_arr, γ_d_arr, γ_l_arr, y_arr, S_arr,
                      istart_arr, istop_arr, N_lines, ::Lorentz, CEF)

KernelAbstractions batched Lorentz line-shape kernel. Each workitem owns one
spectral grid index, loops over all active HITRAN lines, applies each line's
wing-cutoff index window, and accumulates the collision-broadened profile into
`A[I]`.
"""
@kernel function line_shape_batch!(A, @Const(grid), @Const(ν_arr), @Const(γ_d_arr),
                                   @Const(γ_l_arr), @Const(y_arr), @Const(S_arr),
                                   @Const(istart_arr), @Const(istop_arr),
                                   N_lines, ::Lorentz, CEF)
    I = @index(Global, Linear)
    FT = eltype(A)
    acc = zero(FT)
    ν_i = FT(grid[I])
    @inbounds for j in 1:N_lines
        if istart_arr[j] <= I <= istop_arr[j]
            acc += FT(S_arr[j]) * FT(γ_l_arr[j]) /
                   (FT(pi) * (FT(γ_l_arr[j])^2 + (ν_i - FT(ν_arr[j]))^2))
        end
    end
    @inbounds A[I] += acc
end

"""
    line_shape_batch!(A, grid, ν_arr, γ_d_arr, γ_l_arr, y_arr, S_arr,
                      istart_arr, istop_arr, N_lines, ::Voigt, CEF)

KernelAbstractions batched Voigt line-shape kernel. Each workitem owns one
spectral grid index, loops over all active HITRAN lines, applies the line
wing-cutoff window, evaluates the selected complex error function `CEF`, and
accumulates the Voigt profile into `A[I]`.
"""
@kernel function line_shape_batch!(A, @Const(grid), @Const(ν_arr), @Const(γ_d_arr),
                                   @Const(γ_l_arr), @Const(y_arr), @Const(S_arr),
                                   @Const(istart_arr), @Const(istop_arr),
                                   N_lines, ::Voigt, CEF)
    I = @index(Global, Linear)
    FT = eltype(A)
    acc = zero(FT)
    ν_i = FT(grid[I])
    @inbounds for j in 1:N_lines
        if istart_arr[j] <= I <= istop_arr[j]
            acc += FT(S_arr[j]) * FT(cSqrtLn2divSqrtPi) / FT(γ_d_arr[j]) *
                   real(w(CEF, FT(cSqrtLn2) / FT(γ_d_arr[j]) * (ν_i - FT(ν_arr[j])) + im * FT(y_arr[j])))
        end
    end
    @inbounds A[I] += acc
end

#=

Legacy per-line kernels (kept for reference and potential single-line use cases).

=#

"""
    line_shape!(A, grid, ν, γ_d, γ_l, y, S, ::Doppler, CEF)

KernelAbstractions single-line Doppler kernel. Each workitem owns one spectral
grid index and adds this line's Gaussian contribution to `A[I]`.
"""
@kernel function line_shape!(A, @Const(grid), ν, γ_d, γ_l, y, S, ::Doppler, CEF)
    FT = eltype(ν)
    I = @index(Global, Linear)
    @inbounds A[I] += FT(S) * FT(cSqrtLn2divSqrtPi) * exp(-FT(cLn2) * ((FT(grid[I]) - FT(ν)) / FT(γ_d))^2) / FT(γ_d)
end

"""
    line_shape!(A, grid, ν, γ_d, γ_l, y, S, ::Lorentz, CEF)

KernelAbstractions single-line Lorentz kernel. Each workitem owns one spectral
grid index and adds this line's collision-broadened contribution to `A[I]`.
"""
@kernel function line_shape!(A, @Const(grid), ν, γ_d, γ_l, y, S, ::Lorentz, CEF)
    FT = eltype(ν)
    I = @index(Global, Linear)
    @inbounds A[I] += FT(S) * FT(γ_l) / (FT(pi) * (FT(γ_l)^2 + (FT(grid[I]) - FT(ν))^2))
end

"""
    line_shape!(A, grid, ν, γ_d, γ_l, y, S, ::Voigt, CEF)

KernelAbstractions single-line Voigt kernel. Each workitem owns one spectral
grid index, evaluates the selected complex error function `CEF`, and adds the
line's Voigt contribution to `A[I]`.
"""
@kernel function line_shape!(A, @Const(grid), ν, γ_d, γ_l, y, S, ::Voigt, CEF)
    FT = eltype(ν)
    I = @index(Global, Linear)
    @inbounds A[I] += FT(S) * FT(cSqrtLn2divSqrtPi) / FT(γ_d) * real(w(CEF, FT(cSqrtLn2) / FT(γ_d) * (FT(grid[I]) - FT(ν)) + im * FT(y)))
end

"""
    line_shape32!(A, grid, ν, γ_d, γ_l, y, S, ::Doppler, CEF)

Float32-oriented single-line Doppler kernel retained for callers that pass
Float64 host scalars into a Float32 accumulation array. The scalar line
parameters are converted to `eltype(A)` inside the kernel and each workitem
adds this line's Gaussian contribution to one spectral grid point.
"""
@kernel function line_shape32!(A, @Const(grid), ν, γ_d, γ_l, y, S, ::Doppler, CEF)
    FT = eltype(A)
    I = @index(Global, Linear)
    @inbounds A[I] += FT(S) * FT(cSqrtLn2divSqrtPi) * exp(-FT(cLn2) * ((FT(grid[I]) - FT(ν)) / FT(γ_d))^2) / FT(γ_d)
end

"""
    line_shape32!(A, grid, ν, γ_d, γ_l, y, S, ::Lorentz, CEF)

Float32-oriented single-line Lorentz kernel. Each workitem converts scalar
line parameters to `eltype(A)` and adds the collision-broadened line profile
for one spectral grid point.
"""
@kernel function line_shape32!(A, @Const(grid), ν, γ_d, γ_l, y, S, ::Lorentz, CEF)
    FT = eltype(A)
    I = @index(Global, Linear)
    @inbounds A[I] += FT(S) * FT(γ_l) / (FT(pi) * (FT(γ_l)^2 + (FT(grid[I]) - FT(ν))^2))
end

"""
    line_shape32!(A, grid, ν, γ_d, γ_l, y, S, ::Voigt, CEF)

Float32-oriented single-line Voigt kernel. Each workitem evaluates the Voigt
complex-error-function profile in `eltype(A)` and adds the contribution for
one spectral grid point.
"""
@kernel function line_shape32!(A, @Const(grid), ν, γ_d, γ_l, y, S, ::Voigt, CEF)
    FT = eltype(A)
    I = @index(Global, Linear)
    @inbounds A[I] += FT(S) * FT(cSqrtLn2divSqrtPi) / FT(γ_d) * real(w(CEF, FT(cSqrtLn2) / FT(γ_d) * (FT(grid[I]) - FT(ν)) + im * FT(y)))
end

#=

Function to interpolate partition sum for specified isotopologue

=#

"""
    qoft!(M, I, T, T_ref, result)

Compute partition sum ratio Q(T_ref)/Q(T) for temperature corrections.

Uses TIPS 2017 partition function data. Interpolates Q(T) and Q(T_ref) for the given
molecule and isotopologue, then stores Q(T_ref)/Q(T) in result[1].

# Arguments
- `M::Int`: HITRAN molecule ID
- `I::Int`: HITRAN isotopologue ID
- `T`: Target temperature (K)
- `T_ref`: Reference temperature (K)
- `result`: Vector of length 1; result[1] is overwritten with Q(T_ref)/Q(T)

# Throws
- `AssertionError`: If T is outside the TIPS 2017 temperature range for (M, I)
"""
function qoft!(M, I, T, T_ref, result)

    # Get temperature grid
    TT = get_TT(M, I)
    TQ = get_TQ(M, I)

    # Error if out of temperature range
    Tmin = minimum(TT); Tmax = maximum(TT)
    @assert (Tmin < T < Tmax) "TIPS2017: T ($T) must be between $Tmin K and $Tmax K."

    # Interpolate partition sum for specified isotopologue
    interp = DI_CS(TQ, TT)
    Qt = interp(T)
    Qt2 = interp(T_ref)

    # Save the ratio result
    result[1] = Qt2/Qt
end
