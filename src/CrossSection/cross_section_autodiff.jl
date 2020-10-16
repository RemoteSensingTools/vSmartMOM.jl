function test_autodiff()


end

function acs_shorthand(x)
    hitran_data = read_hitran("/home/rjeyaram/RadiativeTransfer/test/helper/CO2.data", mol=2, iso=1, ν_min=6000, ν_max=6400)
    model = make_hitran_model(hitran_data, Doppler(), architecture = Architectures.CPU())
    return absorption_cross_section(model, 6000:0.01:6400, x[1], x[2])
end