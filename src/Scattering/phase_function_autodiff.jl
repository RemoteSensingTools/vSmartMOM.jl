#=
 
This file wraps `compute_aerosol_optical_properties` so that autodiff users and non-autodiff 
users can call the same function with just a keyword argument change. 

=#

""" Function used by auto-differentiation to convert a jacobian result to an AerosolOptics type"""
function convert_jacobian_result_to_aerosol_optics(result)
    
    value = result.value
    derivs = result.derivs[1]

    greek_length = Int64((length(value) - 2)/6)

    α = value[1:greek_length]
    β = value[greek_length + 1 : 2greek_length]
    γ = value[2greek_length + 1 : 3greek_length]
    δ = value[3greek_length + 1 : 4greek_length]
    ϵ = value[4greek_length + 1 : 5greek_length]
    ζ = value[5greek_length + 1 : 6greek_length]

    greek_coefs = GreekCoefs(α, β, γ, δ, ϵ, ζ)

    ω̃ = value[end-1]
    k = value[end]

    return AerosolOptics(greek_coefs=greek_coefs, ω̃=ω̃, k=k, derivs=derivs, fᵗ=eltype(ω̃)(1)) 
end

"""
    $(FUNCTIONNAME)(model::MieModel{FDT})

Reference: Suniti Sanghavi 2014, https://doi.org/10.1016/j.jqsrt.2013.12.015

This function enables user to specify whether to perform auto-differentiation (using either computation type)
Input: MieModel, holding all computation and aerosol properties & autodiff flag (whether to perform 
auto-differentiation)
Output: AerosolOptics, holding all Greek coefficients and Cross-Sectional information
"""
function compute_aerosol_optical_properties(model::MieModel ; autodiff=false)

    # This function takes in the "x-vector" along with the input model so that ForwardDiff will work
    function compute_aerosol_optical_properties_autodiff(x ; model::MieModel = model)

        if length(x) !== 4 
            @error "Must receive four aerosol parameters for auto-differentiation (μ, σ, nᵣ, nᵢ)" x
        end
    
        # Make sure that 𝐱 and model match
        @assert (model.aerosol.size_distribution.μ == log(x[1]))
        @assert (model.aerosol.size_distribution.σ == log(x[2]))
        @assert (model.aerosol.nᵣ == x[3])
        @assert (model.aerosol.nᵢ == x[4])

        # Unpack the model and aerosol 
        @unpack computation_type, aerosol, λ, polarization_type, truncation_type, r_max, nquad_radius, wigner_A, wigner_B = model
        @unpack size_distribution, nᵣ, nᵢ = aerosol

        aerosol_x = Aerosol(LogNormal(log(x[1].value), log(x[2].value)), x[3].value, x[4].value)
        model_x = MieModel(computation_type, aerosol_x, λ, polarization_type, truncation_type, r_max, nquad_radius, wigner_A, wigner_B)
    
        aerosol_optics = compute_aerosol_optical_properties(model_x, ForwardDiff.Dual);
    
        return [aerosol_optics.greek_coefs.α 
                aerosol_optics.greek_coefs.β 
                aerosol_optics.greek_coefs.γ 
                aerosol_optics.greek_coefs.δ 
                aerosol_optics.greek_coefs.ϵ 
                aerosol_optics.greek_coefs.ζ
                aerosol_optics.ω̃
                aerosol_optics.k]
    end

    if (autodiff)

        x = [exp(model.aerosol.size_distribution.μ), 
            exp(model.aerosol.size_distribution.σ), 
            model.aerosol.nᵣ, 
            model.aerosol.nᵢ]

        # Get length of greek coefs
        r, wᵣ = gauleg(model.nquad_radius, 0.0, model.r_max ; norm=true)
        N_max = get_n_max(2 * π * model.r_max/ model.λ)
        greek_length = 2 * N_max - 1

        result = DiffResults.JacobianResult(zeros(6 * greek_length + 2), x)
        ForwardDiff.jacobian!(result, compute_aerosol_optical_properties_autodiff, x);
        return convert_jacobian_result_to_aerosol_optics(result);

    else 
        return compute_aerosol_optical_properties(model)
    end

end