
"Read atmospheric profile (just works for our file, can be generalized"
function read_atmos_profile(file::String, lat::Real, lon::Real, time_idx::Int; g₀=9.8196)

    # Time index must be ∈ [1, 2, 3, 4]
    @assert 1 <= time_idx <= 4 "Time index must be ∈ [1, 2, 3, 4]" 

    # Load in the atmospheric profile
    @timeit "loading file" ds = Dataset(file)

    # See how easy it is to actually extract data? 
    # Note the [:] in the end reads in ALL the data in one step
    file_lats, file_lons = ds["YDim"][:], ds["XDim"][:]
    
    # Convert the input lat/lon to right type
    FT = eltype(file_lats)
    lat, lon = FT(lat), FT(lon)
    
    # Find index (nearest neighbor, one could envision interpolation in space and time!):
    lat_idx, lon_idx = argmin(abs.(file_lats .- lat)), argmin(abs.(file_lons .- lon))

    # Temperature profile
    @timeit "getting T" T = convert(Array{FT,1}, ds["T"][lon_idx, lat_idx,  :, time_idx])

    # Specific humidity profile
    q = convert(Array{FT,1}, ds["QV"][lon_idx, lat_idx, :, time_idx])
    
    # Surface pressure
    psurf = convert(FT, ds["PS"][lon_idx, lat_idx, time_idx])
    
    # AK and BK global attributes (important to calculate pressure half-levels)
    ak, bk = ds.attrib["HDF_GLOBAL.ak"][:], ds.attrib["HDF_GLOBAL.bk"][:]

    # Calculate pressure levels
    p_half = (ak + bk * psurf)
    p_full = (p_half[2:end] + p_half[1:end - 1]) / 2

    # Close the file
    close(ds)

    # Dry and wet mass
    dry_mass = 28.9647e-3  / Nₐ  # in kg/molec, weighted average for N2 and O2
    wet_mass = 18.01528e-3 / Nₐ  # just H2O
    ratio = dry_mass / wet_mass 
    n_layers = length(T)

    # Also get a VMR vector of H2O (volumetric!)
    vmr_h2o = zeros(FT, n_layers, )
    vcd_dry = zeros(FT, n_layers, )
    vcd_h2o = zeros(FT, n_layers, )

    # Now actually compute the layer VCDs
    for i = 1:n_layers 
        Δp = p_half[i + 1] - p_half[i]
        vmr_h2o[i] = q[i] * ratio
        vmr_dry = 1 - vmr_h2o[i]
        M  = vmr_dry * dry_mass + vmr_h2o[i] * wet_mass
        vcd_dry[i] = vmr_dry * Δp / (M * g₀ * 100.0^2)   # includes m2->cm2
        vcd_h2o[i] = vmr_h2o[i] * Δp / (M * g₀ * 100^2)
    end

    # Return the atmospheric profile struct
    return AtmosphericProfile(lat, lon, psurf, T, q, p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o)
end

"Reduce profile dimensions"
function reduce_profile(n::Int, profile::AtmosphericProfile{FT}) where {FT}
    @assert n < length(profile.T)
    @unpack lat, lon, psurf = profile
    # New rough half levels (boundary points)
    a = range(0, maximum(profile.p), length=n + 1)
    # dims = size(σ_matrix)
    # FT = eltype(σ_matrix)
    # σ_matrix_lr = zeros(FT, dims[1], n, dims[3])
    T = zeros(FT, n);
    q = zeros(FT, n);
    p_full = zeros(FT, n);
    p_levels = zeros(FT, n + 1);
    vmr_h2o  = zeros(FT, n);
    vcd_dry  = zeros(FT, n);
    vcd_h2o  = zeros(FT, n);

    for i = 1:n
        ind = findall(a[i] .< profile.p .<= a[i + 1]);
        # σ_matrix_lr[:,i,:] = mean(σ_matrix[:,ind,:], dims=2);
        p_levels[i] = profile.p_levels[ind[1]]
        p_levels[i + 1] = profile.p_levels[ind[end]]
        p_full[i] = mean(profile.p_levels[ind])
        T[i] = mean(profile.T[ind])
        q[i] = mean(profile.q[ind])
        vmr_h2o[i] = mean(profile.vmr_h2o[ind])
        vcd_dry[i] = sum(profile.vcd_dry[ind])
        vcd_h2o[i] = sum(profile.vcd_h2o[ind])
    end

    return AtmosphericProfile(lat, lon, psurf, T, q, p_full, p_levels, vmr_h2o, vcd_dry, vcd_h2o)
end;

# for terrestrial atmospheres 
# psurf in hPa, λ in μm 
function getRayleighLayerOptProp(psurf, λ, depol_fct, vcd_dry) 
    # Total vertical Rayleigh scattering optical thickness 
    tau_scat = 0.00864 * (psurf / 1013.25) * λ^(-3.916 - 0.074 * λ - 0.05 / λ) 
    tau_scat = tau_scat * (6.0 + 3.0 * depol_fct) / (6.0 - 7.0 * depol_fct)
    @show psurf, tau_scat, depol_fct
    Nz = length(vcd_dry)
    τRayl = zeros(Nz)
    k = tau_scat / sum(vcd_dry)
    for i = 1:Nz
        τRayl[i] = k * vcd_dry[i]
    end

    return τRayl
end

# Gaussian distribution on a pressure grid
function getAerosolLayerOptProp(total_τ, p₀, σp, p_half)
    Nz = length(p_half)
    ρ = zeros(Nz)
    for i = 2:Nz
        dp = p_half[i] - p_half[i - 1]
        ρ[i] = (1 / (σp * sqrt(2π))) * exp(-(p_half[i] - p₀)^2 / (2σp^2)) * dp
         # @show ρ[i]  
    end
    Norm = sum(ρ)
    # @show Norm
    τAer  =  (total_τ / Norm) * ρ
    return τAer
end

# computes the composite single scattering parameters (τ, ϖ, Z⁺⁺, Z⁻⁺) for a given atmospheric layer iz for a given Fourier component m
function construct_atm_layer(τRayl, τAer, ϖRayl, ϖAer, fᵗ, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺, τ_abs, arr_type)
    FT = eltype(τRayl)
    # @show FT
    @assert length(τAer) == length(ϖAer) == length(fᵗ) "Sizes don't match"
    
    # @show τRayl , sum(τAer)

    τ = FT(0)
    ϖ = FT(0)
    A = FT(0)
    Z⁺⁺ = similar(Rayl𝐙⁺⁺); 
    Z⁻⁺ = similar(Rayl𝐙⁺⁺);
    
    if (τRayl + sum(τAer)) < eps(FT)
        fill!(Z⁺⁺, 0); fill!(Z⁻⁺, 0);
        return FT(0), FT(1), Z⁺⁺, Z⁻⁺
    end
    # @show τRayl, ϖRayl
    τ += τRayl
    ϖ += τRayl * ϖRayl
    A += τRayl * ϖRayl

    Z⁺⁺ = τRayl * ϖRayl * Rayl𝐙⁺⁺
    Z⁻⁺ = τRayl * ϖRayl * Rayl𝐙⁻⁺

    for i = 1:length(τAer)
        τ   += τAer[i]
        ϖ   += τAer[i] * ϖAer[i]
        # @show τAer[i], ϖAer[i], (1 - fᵗ[i])
        A   += τAer[i] * ϖAer[i] * (1 - fᵗ[i])
        Z⁺⁺ += τAer[i] * ϖAer[i] * (1 - fᵗ[i]) * Aer𝐙⁺⁺[:,:,i]
        Z⁻⁺ += τAer[i] * ϖAer[i] * (1 - fᵗ[i]) * Aer𝐙⁻⁺[:,:,i]
    end
    
    Z⁺⁺ /= A
    Z⁻⁺ /= A
    A /= ϖ
    ϖ /= τ
    
    # Rescaling composite SSPs according to Eqs. A.3 of Sanghavi et al. (2013) or Eqs.(8) of Sanghavi & Stephens (2015)
    # @show A, ϖ
    τ *= (FT(1) - (FT(1) - A) * ϖ)
    ϖ *= ϖ * A / (1 - (1 - A) * ϖ)

    # Adding absorption optical depth / albedo:
    τ_λ = τ_abs .+ τ
    ϖ_λ = (τ .* ϖ) ./ τ_λ
    
    return arr_type(τ_λ), arr_type(ϖ_λ), τ, ϖ, Z⁺⁺, Z⁻⁺  
end

function compute_absorption_profile!(τ_abs::Array{Float64,2}, 
                                     model::AbstractCrossSectionModel,
                                     grid,
                                     profile::AtmosphericProfile,
                                     )

    # pass in the hitran model

    @assert size(τ_abs)[2] == length(profile.p)

    for iz in 1:length(profile.p)

        # Pa -> hPa
        p = profile.p[iz] / 100
        T = profile.T[iz]

        τ_abs[:,iz] = Array(absorption_cross_section(model, grid, p, T)) * profile.vcd_dry[iz] * model.vmr
    end

    return nothing
    
end
