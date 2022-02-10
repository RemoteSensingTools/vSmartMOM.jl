const c₂                 = 1.4387769
const cMassMolIE           = 1.66053873e-30 #kilograms per molecule for unit molec. mass
const cSqrtLn2divSqrtPi  = 0.469718639319144059835 #√(ln2/π)
const cLn2               = 0.6931471805599 #ln2
const cSqrtLn2           = 0.8325546111577 #√(ln2)
const cSqrt2Ln2          = 1.1774100225 #√(2ln2)
const cc_                = 2.99792458e8 #speed of light [m/s]
const cBolts_            = 1.3806503e-23 #Boltzmann const. [J/K]
const p_ref              = 1013.25  # reference pressure [hPa]
const t_ref              = 296.0    # reference temperature [K]
const nm_per_m           = 1.0e7
# Note: ν stands for wavenumber in the following (NOT frequency)
function apply_lineshape_!(Δνᵢ, σᵢ,  # discrete transitions
                λ₀,                 # incident  wavelength [nm]
                # Note:  λ₀ can either be used to denote an individual 
                # wavelength or as a representative (e.g. central) 
                # wavelength for a specific band
                Δν_out,             # Output grid (equidistant)
                σ_out,              # σ at output grid
                pressure::Real,     # actual pressure [hPa]
                temperature::Real,  # Temperature (K)
                molMass; 
                wavelength_flag::Bool=false)
    
    # Notify user of wavelength grid
    if (wavelength_flag)
        @info """
        Note: Rayleigh/Raman Cross-section reported to wavelength grid (nm)
        """  maxlog = 5
    end

    # Max min range (ignoring wing cutoff here)
    grid_max = maximum(Δν_out) 
    grid_min = minimum(Δν_out) 

    # Interpolators from grid bounds to index values
    grid_idx_interp_low  = LinearInterpolation(Δν_out, 1:1:length(Δν_out), extrapolation_bc=1)
    grid_idx_interp_high = LinearInterpolation(Δν_out, 1:1:length(Δν_out), extrapolation_bc=length(Δν_out))

    # Loop through all transition lines:
    for j in eachindex(Δνᵢ)
        # Test that this ν lies within the grid
        if grid_min < Δνᵢ[j] < grid_max

            
            ν = Δνᵢ[j] + nm_per_m/λ₀#13500.0 #Dummy for now #Suniti

            # Compute Doppler HWHM, ν still needs to be supplied, @Suniti?:
            γ_d = ((cSqrt2Ln2 / cc_) * sqrt(cBolts_ / cMassMolIE) * sqrt(temperature) * ν / sqrt(molMass))

            @show γ_d, Δνᵢ[j]

            # line intensity 
            S = σᵢ[j] *  ν^4 #Suniti

            wing_cutoff = 4γ_d 

            # Calculate index range that this transition impacts
            ind_start = Int64(floor(grid_idx_interp_low(Δνᵢ[j] - wing_cutoff)))
            ind_stop  = Int64(ceil(grid_idx_interp_high(Δνᵢ[j] + wing_cutoff)))
            
            # Create views from the result and grid arrays
            result_view   = view(σ_out,  ind_start:ind_stop);
            grid_view     = view(Δν_out, ind_start:ind_stop);

            # Just Doppler broadening for now, can be modified with any line-shape later (using a kernel ideally)
            for I in eachindex(grid_view)
                # If we undersample the line-width, we have to make sure the integral is conserved (almost), TBD
                @inbounds result_view[I] += S * cSqrtLn2divSqrtPi * exp(-cLn2 * (((grid_view[I]) - Δνᵢ[j]) / γ_d)^2) / γ_d
            end
        end
    end
end
