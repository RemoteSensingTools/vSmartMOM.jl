using Revise
using RadiativeTransfer
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM
using Parameters
using LinearAlgebra

# Sets all the "specific" parameters
parameters = vSmartMOM.default_parameters();

# Generates all the derived attributes from above parameters
model = default_model(parameters);

pol_type = model.params.polarization_type;           # Polarization type (IQUV)
obs_geom = model.obs_geom::ObsGeometry;  #, # Solar Zenith, Viewing Zenith, Viewing Azimuthal 
τRayl = model.τRayl;        # Rayleigh optical depth 
    #nAer,                 # Number of aerosol species 
τAer =  model.τAer ;                # Aerosol optical depth and single-scattering albedo
qp_μ =model.qp_μ;
wt_μ =   model.wt_μ  ;      # Quadrature points and weights
Ltrunc =  3  ;            # Trunction length for legendre terms
aerosol_optics = model.aerosol_optics;       # AerosolOptics (greek_coefs, ω̃, k, fᵗ)
GreekRayleigh =  model.greek_rayleigh ;       # Greek coefficients of Rayleigh Phase Function
τ_abs =    model.τ_abs  ;           # nSpec x N
architecture = model.params.architecture;

@unpack obs_alt, sza, vza, vaz = obs_geom   # Observational geometry properties
    
FT = eltype(sza)                    # Get the float-type to use
Nz = length(τRayl)                  # Number of vertical slices
nSpec = size(τ_abs, 1)              # Number of spectral points
μ0 = cosd(sza)                      # μ0 defined as cos(θ); θ = sza
iμ0 = vSmartMOM.nearest_point(qp_μ, μ0)       # Find the closest point to μ0 in qp_μ
arr_type = array_type(architecture)

# Output variables: Reflected and transmitted solar irradiation at TOA and BOA respectively
R = zeros(FT, length(vza), pol_type.n, nSpec)
T = zeros(FT, length(vza), pol_type.n, nSpec)
R_SFI = zeros(FT, length(vza), pol_type.n, nSpec)
T_SFI = zeros(FT, length(vza), pol_type.n, nSpec)

# Copy qp_μ "pol_type.n" times
qp_μN = arr_type(reshape(transpose(repeat(qp_μ, 1, pol_type.n)),pol_type.n*size(qp_μ)[1],1))
#for i = 1:length(qp_μN)
#   @show(i,qp_μN[i]) 
#end
println("Processing on: ", architecture)
println("With FT: ", FT)

#= 
Loop over number of truncation terms =#
SFI = true

m = 1

println("Fourier Moment: ", m)

# Azimuthal weighting
weight = m == 0 ? FT(0.5) : FT(1.0)
# Compute Z-moments of the Rayleigh phase matrix 
# For m>=3, Rayleigh matrices will be 0, can catch with if statement if wanted 
Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, qp_μ, GreekRayleigh, m, arr_type = arr_type);
# Number of aerosols
#@show size(aerosol_optics)
#nBand = length(aerosol_optics)
nAer  = length(aerosol_optics)
# Just for now:
iBand = 1
#nAer, nBand = size(aerosol_optics)
@show nAer#, nBand
dims = size(Rayl𝐙⁺⁺)

# Compute aerosol Z-matrices for all aerosols
Aer𝐙⁺⁺ = arr_type(zeros(FT, (dims[1], dims[2], nAer)))
Aer𝐙⁻⁺ = similar(Aer𝐙⁺⁺)

for i = 1:nAer
    @show aerosol_optics[i,1]
    Aer𝐙⁺⁺[:,:,i], Aer𝐙⁻⁺[:,:,i] = Scattering.compute_Z_moments(pol_type, qp_μ, aerosol_optics[i].greek_coefs, m, arr_type = arr_type)
end

# R and T matrices for Added and Composite Layers for this m
added_layer = vSmartMOM.make_added_layer(FT, arr_type, dims, nSpec) 
composite_layer = vSmartMOM.make_composite_layer(FT, arr_type, dims, nSpec)
I_static = Diagonal(arr_type(Diagonal{FT}(ones(dims[1]))));
scattering_interface = vSmartMOM.ScatteringInterface_00()
τ_sum = zeros(nSpec) #Suniti: declaring τ_sum to be of length nSpec
τ_λ = zeros(nSpec)
iz = Nz
if iz==1
    τ_sum = τ_λ
else
    τ_sum = τ_sum + τ_λ     
end
#@show(iz, Nz)
# Construct the atmospheric layer
# From Rayleigh and aerosol τ, ϖ, compute overall layer τ, ϖ
#@timeit "Constructing" 

τ_λ, ϖ_λ, τ, ϖ, Z⁺⁺, Z⁻⁺ = vSmartMOM.construct_atm_layer(τRayl[iz], τAer[:,iz], aerosol_optics, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, Aer𝐙⁺⁺, Aer𝐙⁻⁺, τ_abs[:,iz], arr_type)
#@show(τ_λ)
#@show(ϖ_λ)
#@show(τ)
#@show(ϖ)
#sleep(5)
#for i=1:size(Z⁺⁺)[1]
#    @show(i,Z⁺⁺[i,:])
#end
# τ * ϖ should remain constant even though they individually change over wavelength
# @assert all(i -> (i ≈ τ * ϖ), τ_λ .* ϖ_λ)
# Compute doubling number
dτ_max = minimum([τ * ϖ, FT(0.01) * minimum(qp_μ)])
dτ, ndoubl = vSmartMOM.doubling_number(dτ_max, τ * ϖ) #Suniti
#@show(ndoubl, dτ_max, τ)
# Compute dτ vector
dτ_λ = arr_type(τ_λ ./ (FT(2)^ndoubl))
expk = exp.(-dτ_λ /qp_μ[iμ0]) #Suniti
#@show(τ_λ, dτ_λ.*FT(2)^ndoubl)
#@show(τ, dτ*FT(2)^ndoubl,dτ, dτ_λ )
#@show(expk, exp.(-dτ /qp_μ[iμ0]))
#@show τ_sum
#@show dτ_λ, dτ
#scatter = true
vSmartMOM.elemental!(pol_type, SFI, iμ0, τ_sum, dτ_λ, dτ, ϖ_λ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, true, qp_μ, wt_μ, added_layer,  I_static, arr_type, architecture)
vSmartMOM.doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
#added_layer_DNI = vSmartMOM.make_added_layer(FT, arr_type, dims, nSpec) 
#composite_layer_DNI = vSmartMOM.make_composite_layer(FT, arr_type, dims, nSpec)
#vSmartMOM.elemental!(pol_type, false, iμ0, τ_sum, dτ_λ, dτ, ϖ_λ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, true, qp_μ, wt_μ, added_layer_DNI,  I_static, arr_type, architecture)