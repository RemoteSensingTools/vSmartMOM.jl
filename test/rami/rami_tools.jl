function add_aerosols!(rami_atmosphere, params)# TODO We still need to check the two fractions of coarse and fine! Right now, only one is
    if length(rami_atmosphere["aerosols"]) > 0# 
        if startswith(rami_atmosphere["aerosols"][1]["name"], "D")
            @show "Desert Aerosol"
            μ_fine   = 0.0478666
            σ_fine   = 1.87411
            μ_coarse = 0.604127
            σ_coarse = 1.75172
            # Scaled so that sum is 1
            n_coarse = 0.0033219635
            
            # These refractive indices depend on the band. 
            # TODO: We need to pull those from the table actually, as they are band dependent 
            # (in vSmartMOM, we so far only have one ni/nr per type).
            # See https://rami-benchmark.jrc.ec.europa.eu/_www/RAMI4ATM/down/RAMI4ATM_aerosols_v1.0/refractive_index/continental.txt
            # Need to generalize
            
            # TODO: Make it wavelength dependent (and reference scenario fixed)
            nᵣ = 1.477538814814815
            nᵢ = 0.004342592592592592

        elseif startswith(rami_atmosphere["aerosols"][1]["name"], "C")
            @show "Continental Aerosol"
            #TODO: Change, not just one refractice index
            μ_fine   = 0.0807989
            σ_fine   = 1.50180
            μ_coarse = 0.682651
            σ_coarse = 2.10400
            # Scaled so that sum is 1
            n_coarse   = 0.00046374026257

            # TODO: Make it wavelength dependent (and reference scenario fixed)
            nᵣ = 1.4434925925925925
            nᵢ = 0.0015797

        else
            println("weird aerosol here")
        end

        # TODO: Suniti, please double-check

        # size_distribution = LogNormal(log(μ), log(σ))
        size_distribution = MixtureModel(LogNormal,
                            [(log(μ_fine), log(σ_fine)), (log(μ_coarse), log(σ_coarse))],
                            [1-n_coarse, n_coarse])

        # Create the aerosol(s)
        aero = vSmartMOM.CoreRT.Aerosol(size_distribution, nᵣ, nᵢ)

        # Need to do the pressure peak more carefully
        # TODO: We need to just put it into 1-2 layers...
        RT_aerosol = vSmartMOM.CoreRT.RT_Aerosol(aero, rami_atmosphere["aerosols"][1]["tau_550"], Uniform(795.0,1013.0))

        # Assemble scattering parameters
        scattering_params = vSmartMOM.CoreRT.ScatteringParameters([RT_aerosol], 30.0, 2500, 0.550, vSmartMOM.Scattering.NAI2())

        params.scattering_params = scattering_params
        @show params.scattering_params.rt_aerosols[1].τ_ref
    end
end

# Get q from vmr:
q = -1/((dry_mass - dry_mass./vmr_h2o)/wet_mass - 1)