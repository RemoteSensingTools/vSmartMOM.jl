const c₂                 = 1.4387769
const cMassMol           = 1.66053873e-27 #grams per molecule for unit molec. mass
const cSqrtLn2divSqrtPi  = 0.469718639319144059835 #√(ln2/π)
const cLn2               = 0.6931471805599 #ln2
const cSqrtLn2           = 0.8325546111577 #√(ln2)
const cSqrt2Ln2          = 1.1774100225 #√(2ln2)
const cc_                = 2.99792458e8 #speed of light [m/s]
const cBolts_            = 1.3806503e-23 #Boltzmann const. [J/K]
const p_ref              = 1013.25  # reference pressure [hPa]
const t_ref              = 296.0    # reference temperature [K]
const nm_per_m           = 1.0e7

function get_n₀_n₁(ieJ₁⁺,Δ)
    n₁_ = 1:size(ieJ₁⁺,3);
    n₀_  = n₁_ .+ Δ  ;
    # Find valid indices:
    sub = findall(1 .≤ n₀_ .≤ size(ieJ₁⁺,3));
    n₁ = n₁_[sub[1]]:n₁_[sub[end]]
    n₀ = n₀_[sub[1]]:n₀_[sub[end]]
    return n₀, n₁
end
# Currently assuming same T for all vertical atmospheric layers (so that a uniform Raman wavelength grid can be assumed for rt_interactions)
function getRamanAtmoConstants(ν̃::FT, T::FT) where FT
    n2 = InelasticScattering.getMolecularConstants(InelasticScattering.N₂(), (0.8));
    compute_effective_coefficents!(ν̃, T, n2)
    compute_energy_levels!(n2)
    compute_σ_Rayl_coeff!(n2)
    compute_σ_Rayl_VibRaman_coeff_hires!(T, n2)
    compute_σ_VibRaman_coeff!(T, n2)
    compute_σ_RoVibRaman_coeff!(T, n2)

    o2 = InelasticScattering.getMolecularConstants(InelasticScattering.O₂(), (0.2));
    compute_effective_coefficents!(ν̃, T, o2)
    compute_energy_levels!(o2)
    compute_σ_Rayl_coeff!(o2)
    compute_σ_Rayl_VibRaman_coeff_hires!(T, o2)
    compute_σ_VibRaman_coeff!(T, o2)
    compute_σ_RoVibRaman_coeff!(T, o2)

    return n2,o2
end

function compute_ϖ_Cabannes(RS_type::noRS, λ₀)
    RS_type.ϖ_Cabannes = 1.0;
    return RS_type.ϖ_Cabannes;
end

function compute_ϖ_Cabannes(RS_type::Union{RRS, VS_0to1, VS_1to0, RRS_plus, VS_0to1_plus, VS_1to0_plus}, λ₀, n2, o2)
    ν₀ = 1e7/λ₀;

    σ_elastic =  n2.vmr * n2.effCoeff.σ_Rayl_coeff + o2.vmr * o2.effCoeff.σ_Rayl_coeff 
    σ_elastic *= ν₀^4
    
    σ_RRS =  n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4)' * n2.effCoeff.σ_RoRaman_coeff_JtoJp2
    σ_RRS += n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4)' * n2.effCoeff.σ_RoRaman_coeff_JtoJm2

    σ_RRS += o2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4)' * n2.effCoeff.σ_RoRaman_coeff_JtoJp2
    σ_RRS += o2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4)' * n2.effCoeff.σ_RoRaman_coeff_JtoJm2

    σ_RVRS =  n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2).^4)' * n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2
    σ_RVRS += n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2).^4)' * n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2
    σ_RVRS += n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2).^4)' * n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2
    σ_RVRS += n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2).^4)' * n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2
    σ_RVRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2).^4)' * o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2
    σ_RVRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2).^4)' * o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2
    σ_RVRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2).^4)' * o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2
    σ_RVRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2).^4)' * o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2

    σ_VRS =  n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires).^4)' * n2.effCoeff.σ_VibRaman_coeff_0to1_hires
    σ_VRS += n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires).^4)' * n2.effCoeff.σ_VibRaman_coeff_1to0_hires
    σ_VRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires).^4)' * o2.effCoeff.σ_VibRaman_coeff_0to1_hires
    σ_VRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires).^4)' * o2.effCoeff.σ_VibRaman_coeff_1to0_hires    

    #RS_type.ϖ_Cabannes = σ_elastic/(σ_VRS+σ_RVRS+σ_RRS+σ_elastic);
    ϖ_Cabannes = σ_elastic/(σ_RRS+σ_elastic);
    return ϖ_Cabannes;
end

# Note: ν stands for wavenumber in the following (NOT frequency)
function apply_lineshape!(Δνᵢ, σᵢ,  # discrete transitions
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
    σ_out .= 0;
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

    S_sum=0.0
    # Loop through all transition lines:
    for j in eachindex(Δνᵢ)
        # Test that this ν lies within the grid
        if grid_min < Δνᵢ[j] < grid_max
            
            ν = Δνᵢ[j] + nm_per_m/λ₀#13500.0 #Dummy for now #Suniti

            # Compute Doppler HWHM, ν still needs to be supplied, @Suniti?:
            γ_d = ((cSqrt2Ln2 / cc_) * sqrt(cBolts_ / cMassMol) * sqrt(temperature) * ν / sqrt(molMass))

            # line intensity 
            S = σᵢ[j] *  ν^4 #Suniti
            S_sum += S
            #@show γ_d, Δνᵢ[j], S, ν, λ₀, cMassMol

            wing_cutoff = 2γ_d 

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
                #@show grid_view[I], result_view[I]
            end
        end
    end
    #dν = Δν_out[2]-Δν_out[1] 
    nothing  
    #@show S_sum, sum(σ_out)*dν
end

function apply_gridlines!(Δνᵢ, σᵢ,  # discrete transitions
    λ₀,                 # incident  wavelength [nm]
    # Note:  λ₀ can either be used to denote an individual 
    # wavelength or as a representative (e.g. central) 
    # wavelength for a specific band
    ν_in,             # Output grid (equidistant)
    σ_out;              # σ at output grid
    
    #pressure::Real,     # actual pressure [hPa]
    #temperature::Real,  # Temperature (K)
    #molMass; 
    wavelength_flag::Bool=false)    
    σ_out .= 0;
    # Notify user of wavelength grid
    if (wavelength_flag)
        @info """
        Note: Rayleigh/Raman Cross-section reported to wavelength grid (nm)
        """  maxlog = 5
    end
    Δν_in = ν_in .- nm_per_m/λ₀

    # Max min range (ignoring wing cutoff here)
    grid_max = maximum(Δν_in) 
    grid_min = minimum(Δν_in) 

    # Interpolators from grid bounds to index values
    #grid_idx_interp_low  = LinearInterpolation(Δν_out, 1:1:length(Δν_out), extrapolation_bc=1) 
    #grid_idx_interp_high = LinearInterpolation(Δν_out, 1:1:length(Δν_out), extrapolation_bc=length(Δν_out))
    S_sum = 0.0
    # Loop through all transition lines:
    for j in eachindex(Δνᵢ)
        #@show(grid_min, Δνᵢ[j], grid_max)
        # Test that this ν lies within the grid
        if grid_min < Δνᵢ[j] < grid_max
            
            ν = Δνᵢ[j] + nm_per_m/λ₀#13500.0 #Dummy for now #Suniti

            # Compute Doppler HWHM, ν still needs to be supplied, @Suniti?:
            #γ_d = ((cSqrt2Ln2 / cc_) * sqrt(cBolts_ / cMassMol) * sqrt(temperature) * ν / sqrt(molMass))
            # line intensity 
            S = σᵢ[j] *  ν^4 #Suniti
            #@show γ_d, Δνᵢ[j], S, ν, λ₀, cMassMol
            S_sum += S
            #wing_cutoff = 2γ_d 
            i=argmin(abs.(Δνᵢ[j].-Δν_in))
            #@show i, Δνᵢ[j]-Δν_in[i]    
            if Δν_in[i]<Δνᵢ[j]
                ind_start = i
                ind_stop  = i+1
            else
                ind_start = i-1
                ind_stop  = i
            end         
            # Calculate index range that this transition impacts
            #ind_start = Int64(floor(grid_idx_interp_low(Δνᵢ[j] - wing_cutoff)))
            #ind_stop  = Int64(ceil(grid_idx_interp_high(Δνᵢ[j] + wing_cutoff)))

            # Create views from the result and grid arrays
            result_view   = view(σ_out,  ind_start:ind_stop);
            grid_view     = view(Δν_in, ind_start:ind_stop);

            # Just Doppler broadening for now, can be modified with any line-shape later (using a kernel ideally)
            for I in eachindex(grid_view)
                # If we undersample the line-width, we have to make sure the integral is conserved (almost), TBD
                @inbounds result_view[I] += S/2.# * cSqrtLn2divSqrtPi * exp(-cLn2 * (((grid_view[I]) - Δνᵢ[j]) / γ_d)^2) / γ_d
                #@show grid_view[I], result_view[I]
            end
        end
    end
    #@show S_sum, sum(σ_out)
    nothing
end


#function compute_optical_Rayl!(grid_out,atmo_σ_Rayl, λ₀, n2, o2)
function compute_optical_Rayl(λ₀, n2, o2)
    atmo_σ_Rayl = n2.vmr * n2.effCoeff.σ_Rayl_coeff #σ_out #cross section in cm^2
    atmo_σ_Rayl += o2.vmr * o2.effCoeff.σ_Rayl_coeff #cross section in cm^2
    atmo_σ_Rayl *= (nm_per_m/λ₀)^4
    #plot(1.e7/λ₀ .+ grid_out,atmo_σ_Rayl*1.e40)
    return atmo_σ_Rayl;
end

function compute_optical_RS!(RS_type::Union{RRS, RRS_plus}, grid_in, λ₀, n2, o2)
    #plotly()
    # grid_in is a uniform wavenumber grid covering the entire band spectrum 
    # TMP: grid_in = nm_per_m/λ₀.+collect(-250:0.002:250) #this is a wavenumber grid
    # get_greek_raman!(RS_type, n2, o2)
    
    σ_out = similar(grid_in);
    atmo_σ_RRS_JtoJp2 = similar(grid_in);
    atmo_σ_RRS_JtoJm2 = similar(grid_in);     
    σ_tmp = similar(grid_in);

    # N2
    #apply_lineshape!(n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2, n2.effCoeff.σ_RoRaman_coeff_JtoJp2,  λ₀, collect(grid_out), σ_out, 1, 300.0, 28);
    apply_gridlines!(n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2, n2.effCoeff.σ_RoRaman_coeff_JtoJp2,  λ₀, collect(grid_in), σ_out);
    atmo_σ_RRS_JtoJp2 = n2.vmr * σ_out #cross section in cm^2
    #@show length(atmo_σ_RRS_JtoJp2[atmo_σ_RRS_JtoJp2.>0])
    #for I in eachindex(grid_out)
    #    @show grid_out[I], σ_out[I]
    #end

    #apply_lineshape!(n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2, n2.effCoeff.σ_RoRaman_coeff_JtoJm2, λ₀, collect(grid_out), σ_out, 1, 300.0, 40);
    apply_gridlines!(n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2, n2.effCoeff.σ_RoRaman_coeff_JtoJm2, λ₀, collect(grid_in), σ_out);
    atmo_σ_RRS_JtoJm2 = n2.vmr * σ_out #cross section in cm^2
    #@show length(atmo_σ_RRS_JtoJm2[atmo_σ_RRS_JtoJm2.>0])
    # O2
    #apply_lineshape!(o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2, o2.effCoeff.σ_RoRaman_coeff_JtoJp2, λ₀, collect(grid_out), σ_out, 1, 300.0, 28);
    apply_gridlines!(o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2, o2.effCoeff.σ_RoRaman_coeff_JtoJp2, λ₀, grid_in, σ_out);
    atmo_σ_RRS_JtoJp2 += o2.vmr * σ_out #cross section in cm^2
    #@show length(atmo_σ_RRS_JtoJp2[atmo_σ_RRS_JtoJp2.>0])

    #apply_lineshape!(o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2, o2.effCoeff.σ_RoRaman_coeff_JtoJm2, λ₀, collect(grid_out), σ_out, 1, 300.0, 40);
    apply_gridlines!(o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2, o2.effCoeff.σ_RoRaman_coeff_JtoJm2, λ₀, grid_in, σ_out);
    atmo_σ_RRS_JtoJm2 += o2.vmr * σ_out #cross section in cm^2
    #@show length(atmo_σ_RRS_JtoJm2[atmo_σ_RRS_JtoJm2.>0])

    σ_tmp .= atmo_σ_RRS_JtoJm2 .+ atmo_σ_RRS_JtoJp2
    atmo_σ_RRS = σ_tmp[σ_tmp.>0]
    #finding all indices of σ_out (and hence of ν_in) that have finite (non-zero) values
    index_ramangrid_out = findall(x->x in σ_tmp[σ_tmp.>0],σ_tmp)
    if (nm_per_m/λ₀>grid_in[1] && nm_per_m/λ₀<grid_in[end])
        index_ramangrid_out .-= argmin(abs.(grid_in .- nm_per_m/λ₀))
    end 
    #for I in eachindex(atmo_σ_RRS)
    #    @show grid_in[argmin(abs.(grid_in .- nm_per_m/λ₀))+index_ramangrid_out[I]], index_ramangrid_out[I], atmo_σ_RRS[I]
    #end    
    return index_ramangrid_out, atmo_σ_RRS;
    #plot(grid_out, atmo_σ_RRS_JtoJp2, yscale=:log10)
    #plot(1.e7/λ₀ .+ grid_out, atmo_σ_RRS_plot*1.e40)
end

function compute_optical_RS!(RS_type::Union{VS_0to1, VS_0to1_plus}, grid_in, λ₀, n2, o2)
    #plotly()
    #get_greek_raman(RS_type, n2, o2)
    #compute_ϖ_Cabannes!(RS_type, λ₀, n2, o2)

    #@show n2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires[0], o2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires[0]
    #νᵣ = 0.5*(n2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires[0] + o2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires[0])
    
    # TMP: grid_in = nm_per_m/λ₀ .+ collect((νᵣ-750):0.002:(νᵣ+750))
    σ_out = similar(grid_in);
    #atmo_σ_VRS_0to1 = similar(grid_in);
    #atmo_σ_RVRS_0to1 = similar(grid_in);
    σ_tmpVRS = similar(grid_in);
    σ_tmpRVRS = similar(grid_in);
    # N2
    xin = [n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2.parent; n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2.parent]
    yin = [n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2.parent; n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2.parent];
    #apply_lineshape!(xin, yin, λ₀, collect(grid_out), σ_out, 1, 300.0, 40);
    apply_gridlines!(xin[abs.(xin).>0], yin[abs.(xin).>0], λ₀, grid_in, σ_out);
    σ_tmpRVRS = n2.vmr * σ_out #cross section in cm^2
    #for i in 1:length(xin)
    #    @show i, xin[i], yin[i]
    #end
    xin = n2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires
    yin = n2.effCoeff.σ_VibRaman_coeff_0to1_hires
    #for i in 0:length(xin)-1
    #    @show i, xin[i], yin[i]
    #end
    #apply_lineshape!(xin, yin, λ₀, collect(grid_out), σ_out, 1, 300.0, 40);
    apply_gridlines!(xin, yin, λ₀, grid_in, σ_out);
    σ_tmpVRS = n2.vmr * σ_out #cross section in cm^2

    # O2
    xin = [o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2.parent; o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2.parent]
    yin = [o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2.parent; o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2.parent];
    #apply_lineshape!(xin, yin, λ₀, collect(grid_out), σ_out, 1, 300.0, 40);
    apply_gridlines!(xin[abs.(xin).>0], yin[abs.(xin).>0], λ₀, grid_in, σ_out);
    σ_tmpRVRS += o2.vmr * σ_out #cross section in cm^2
    #for i in 1:length(σ_out)
    #    @show i, 1.e7/λ₀ .+ grid_out[i], σ_out[i]
    #end
    xin = o2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires
    yin = o2.effCoeff.σ_VibRaman_coeff_0to1_hires
    #apply_lineshape!(xin, yin, λ₀, collect(grid_out), σ_out, 1, 300.0, 40);
    apply_gridlines!(xin, yin, λ₀, grid_in, σ_out);
    σ_tmpVRS += o2.vmr * σ_out #cross section in cm^2

    atmo_σ_VRS_0to1 = σ_tmpVRS[σ_tmpVRS.>0]
    #finding all indices of σ_out (and hence of ν_in) that have finite (non-zero) values
    index_VRSgrid_out = findall(x->x in σ_tmpVRS[σ_tmpVRS.>0],σ_tmpVRS)

    atmo_σ_RVRS_0to1 = σ_tmpRVRS[σ_tmpRVRS.>0]
    #finding all indices of σ_out (and hence of ν_in) that have finite (non-zero) values
    index_RVRSgrid_out = findall(x->x in σ_tmpRVRS[σ_tmpRVRS.>0],σ_tmpRVRS)
    #plot(grid_out, atmo_σ_RRS_JtoJp2, yscale=:log10)
    #plot(grid_in, σ_tmpRVRS*1.e40)
    #plot!(grid_in, σ_tmpVRS*1.e40)
    return index_VRSgrid_out, atmo_σ_VRS_0to1, index_RVRSgrid_out, atmo_σ_RVRS_0to1;
end


function compute_optical_RS!(RS_type::Union{VS_1to0, VS_1to0_plus}, grid_in, λ₀, n2, o2)
    #plotly()
    get_greek_raman(RS_type, n2, o2)
    compute_ϖ_Cabannes!(RS_type, λ₀, n2, o2)
    @show n2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires[0], o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires[0]
    νᵣ = 0.5*(n2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires[0] + o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires[0])
        
    # TMP: grid_in = nm_per_m/λ₀ + collect((νᵣ-750):0.002:(νᵣ+750))
    σ_out = similar(collect(grid_in));        
    atmo_σ_VRS_1to0 = similar(grid_in);
    atmo_σ_RVRS_1to0 = similar(grid_in);
    σ_tmpVRS = similar(grid_in);
    σ_tmpRVRS = similar(grid_in);

    @show n2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires[1], o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires[1]
    #νᵣ = 0.5*(n2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires[1] + o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires[1])
    #grid_out = (νᵣ-750):0.002:(νᵣ+750)
    #σ_out = similar(collect(grid_out));
    # N2
    xin = [n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2.parent; n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2.parent]
    yin = [n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2.parent; n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2.parent];
    #apply_lineshape!(xin, yin, λ₀, collect(grid_out), σ_out, 1, 300.0, 40);
    apply_gridlines!(xin[abs.(xin).>0], yin[abs.(xin).>0], λ₀, grid_in, σ_out);
    σ_RVRStmp= n2.vmr * σ_out #cross section in cm^2

    xin = n2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires
    yin = n2.effCoeff.σ_VibRaman_coeff_1to0_hires
    #apply_lineshape!(xin, yin, λ₀, collect(grid_out), σ_out, 1, 300.0, 40);
    apply_gridlines!(xin, yin, λ₀, grid_in, σ_out);
    σ_VRStmp = n2.vmr * σ_out #cross section in cm^2

    # O2
    xin = [o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2.parent; o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2.parent]
    yin = [o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2.parent; o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2.parent];
    #apply_lineshape!(xin, yin, λ₀, collect(grid_out), σ_out, 1, 300.0, 40);
    
    apply_gridlines!(xin[abs.(xin).>0], yin[abs.(xin).>0], λ₀, grid_in, σ_out);
    σ_RVRStmp += o2.vmr * σ_out #cross section in cm^2
    
    xin = o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires
    yin = o2.effCoeff.σ_VibRaman_coeff_1to0_hires
    #apply_lineshape!(xin, yin, λ₀, collect(grid_out), σ_out, 1, 300.0, 40);
    apply_gridlines!(xin, yin, λ₀, grid_in, σ_out);
    σ_VRStmp += o2.vmr * σ_out #cross section in cm^2

    atmo_σ_VRS_1to0 .= σ_tmpVRS(σ_VRStmp.>0)
    #finding all indices of σ_out (and hence of ν_in) that have finite (non-zero) values
    index_VRSgrid_out = findall(x->x in σ_tmpVRS[σ_tmpVRS.>0],σ_tmpVRS)

    atmo_σ_RVRS_1to0 .= σ_tmpRVRS(σ_RVRStmp.>0)
    #finding all indices of σ_out (and hence of ν_in) that have finite (non-zero) values
    index_RVRSgrid_out = findall(x->x in σ_tmpRVRS[σ_tmpRVRS.>0],σ_tmpRVRS)

    return index_VRSgrid_out, atmo_σ_VRS_1to0, index_RVRSgrid_out, atmo_σ_RVRS_1to0;
    #plot(grid_out, atmo_σ_RRS_JtoJp2, yscale=:log10)
    #plot(grid_in, σ_RVRStmp*1.e40)
    #plot!(grid_in, σ_VRStmp*1.e40)
end

"""
    $(FUNCTIONNAME)(depol)
Returns the greek coefficients (as [`GreekCoefs`](@ref)) of the Rayleigh phase function given 
depolarization value. 
- `depol` Depolarization (best use 0 as default )
"""
function get_greek_raman!(RS_type::noRS, n2, o2)
    return nothing
end

# the following applies to both rovibrational and rotational Raman scattering (by both N2 and O2)
function get_greek_raman(RS_type::Union{RRS, RRS_plus, VS_0to1, VS_0to1_plus, VS_1to0, VS_1to0_plus}, 
                            n2, o2)
    depol = n2.effCoeff.rho_depol_RotRaman
    FT = eltype(depol)

    # Rayleigh Greek Parameters
    dpl_p = (1 - depol)  / (1 + depol / 2)
    #dpl_q = (1 + depol)  / (1 - depol)
    dpl_r = (1 - 2depol) / (1 - depol)
  
    α  =  FT[0.0, 0.0,             3dpl_p]
    β  =  FT[1.0, 0.0,             0.5 * dpl_p]
    γ  =  FT[0.0, 0.0,             dpl_p * sqrt(1.5)] 
    δ  =  FT[0.0, dpl_p * dpl_r * 1.5, 0.0] 
    ϵ  =  FT[0.0, 0.0,             0.0] 
    ζ  =  FT[0.0, 0.0,             0.0]
    return GreekCoefs(α, β, γ, δ, ϵ, ζ);
    #return nothing
end

function get_greek_raman_VS(RS_type::Union{VS_0to1, VS_0to1_plus, VS_1to0, VS_1to0_plus}, 
                            in_molec)
    
    depol = in_molec.effCoeff.rho_depol_VibRaman
    
    FT = eltype(depol)
    # Rayleigh Greek Parameters
    dpl_p = (1 - depol)  / (1 + depol / 2)
    #dpl_q = (1 + depol)  / (1 - depol)
    dpl_r = (1 - 2depol) / (1 - depol)
  
    α  =  FT[0.0, 0.0,             3dpl_p]
    β  =  FT[1.0, 0.0,             0.5 * dpl_p]
    γ  =  FT[0.0, 0.0,             dpl_p * sqrt(1.5)] 
    δ  =  FT[0.0, dpl_p * dpl_r * 1.5, 0.0] 
    ϵ  =  FT[0.0, 0.0,             0.0] 
    ζ  =  FT[0.0, 0.0,             0.0]
    return GreekCoefs(α, β, γ, δ, ϵ, ζ);
    #return nothing
end

function compute_Rayl_depol(n2, o2)
    depol = (n2.vmr*n2.effCoeff.rho_depol_Rayl + o2.vmr*o2.effCoeff.rho_depol_Rayl)/(n2.vmr+o2.vmr)
    return depol
end


function computeRamanZλ!(RS_type::Union{RRS_plus,RRS}, pol_type, qp_μ, m, arr_type)
    RS_type.Z⁺⁺_λ₁λ₀, RS_type.Z⁻⁺_λ₁λ₀ =  Scattering.compute_Z_moments(pol_type, 
                                        qp_μ, 
                                        RS_type.greek_raman, 
                                        m, 
                                        arr_type = arr_type);
    nothing
end

function computeRamanZλ!(RS_type::Union{noRS_plus, noRS}, pol_type, qp_μ, m, arr_type)
    nothing
end

function computeRamanZλ!(RS_type::AbstractRamanType, pol_type, qp_μ, m, arr_type)
    RS_type.Z⁺⁺_λ₁λ₀, RS_type.Z⁻⁺_λ₁λ₀ = Scattering.compute_Z_moments(pol_type, 
                                        qp_μ, 
                                        RS_type.greek_raman, 
                                        m, 
                                        arr_type = arr_type);
    RS_type.Z⁺⁺_λ₁λ₀_VS_n2, RS_type.Z⁻⁺_λ₁λ₀_VS_n2 = 
                    Scattering.compute_Z_moments(pol_type, 
                                            Array(qp_μ), 
                                            RS_type.greek_raman_VS_n2, 
                                            m, 
                                            arr_type = arr_type);
    RS_type.Z⁺⁺_λ₁λ₀_VS_o2, RS_type.Z⁻⁺_λ₁λ₀_VS_o2 = 
                    Scattering.compute_Z_moments(pol_type, 
                                        Array(qp_μ), 
                                        RS_type.greek_raman_VS_o2, 
                                        m, 
                                        arr_type = arr_type);      
    nothing
end




