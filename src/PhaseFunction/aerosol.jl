# This file currently contains just an initializer for UnivariateAerosol but it can also 
# contain helper functions in the future (combining univariate distributions, doing more 
# complex tasks with aerosol distributions, etc.)

""" 
    $(FUNCTIONNAME)(size_distribution::ContinuousUnivariateDistribution, r_max, nquad_radius::Int, nᵣ, nᵢ)

Convenience function to create a Univariate Aerosol, given the size parameters 
"""
function make_univariate_aerosol(size_distribution::ContinuousUnivariateDistribution, 
                                 r_max, nquad_radius::Int, 
                                 nᵣ, nᵢ)

    return UnivariateAerosol(size_distribution, r_max, nquad_radius, nᵣ, nᵢ)
end

