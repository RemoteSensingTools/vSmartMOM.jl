using InelasticScattering
using Plots
λ = 440*1.e-7 #cm
ν̃ = 1/λ
T = 300 #k

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

plotLegend = false
plot_coefficients = 1
plot_cross_sections = 0
if plot_coefficients==1
    @show "YEAH"
    # N2
    # Rayleigh scattering, pure vibrational Raman scattering
    x = [n2.effCoeff.Δν̃_VibRaman_coeff_0to1 n2.effCoeff.Δν̃_VibRaman_coeff_1to0 0];
    y = [n2.effCoeff.σ_VibRaman_coeff_0to1 n2.effCoeff.σ_VibRaman_coeff_1to0 n2.effCoeff.σ_Rayl_coeff];
    @show(ν̃)
    @show(x)
    @show(y)
    #plot(x, y*1e50, color=:blue, seriestype=:sticks, legend=plotLegend)

    # Rayleigh hires scattering
    #x = [n2.effCoeff.Δν̃_Rayl_coeff_hires.-n2.effCoeff.Δν̃_Rayl_coeff_hires[0]];
    #y = [n2.effCoeff.σ_Rayl_coeff_hires];
    #@show(ν̃)
    #@show(x)
    #@show(y)
    #plot(x, y*1e50, color=:blue, seriestype=:sticks, legend=plotLegend)

    # Rotational Raman scattering
    x = [n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2 n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2];
    y = [n2.effCoeff.σ_RoRaman_coeff_JtoJp2 n2.effCoeff.σ_RoRaman_coeff_JtoJm2];
    @show(x)
    @show(y)
    #plot!(x, y * 1e50, color=:blue, seriestype=:sticks, legend=plotLegend)

    #Hires vibrational scattering
    x = [n2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires]# n2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires];
    y = [n2.effCoeff.σ_VibRaman_coeff_0to1_hires]# n2.effCoeff.σ_VibRaman_coeff_1to0_hires];
    @show(x)
    @show(y)
    plot(x, y*1e50, color=:blue, seriestype=:sticks, legend=plotLegend)

    # Rovibrational scattering (P)
    x = [n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2 n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2];
    y = [n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2 n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2];
    @show(x)
    @show(y)
    #plot!(x, y*1e50, color=:blue, seriestype=:sticks, legend=plotLegend)

    # Rovibrational scattering (R)
    x = [n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2 n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2];
    y = [n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2 n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2];
    @show(x)
    @show(y)
    plot!(x, y*1e50, color=:blue, seriestype=:sticks, legend=plotLegend)

    # O2
    # Rayleigh scattering, pure vibrational Raman scattering
    x = [o2.effCoeff.Δν̃_VibRaman_coeff_0to1 o2.effCoeff.Δν̃_VibRaman_coeff_1to0 0];
    y = [o2.effCoeff.σ_VibRaman_coeff_0to1 o2.effCoeff.σ_VibRaman_coeff_1to0 o2.effCoeff.σ_Rayl_coeff];
    @show(ν̃)
    @show(x)
    @show(y)
    #plot!(x, y*1e50, color=:orange, seriestype=:sticks, legend=plotLegend)

    # Rotational Raman scattering
    x = [o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2 o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2];
    y = [o2.effCoeff.σ_RoRaman_coeff_JtoJp2 o2.effCoeff.σ_RoRaman_coeff_JtoJm2];
    @show(x)
    @show(y)
    #plot!(x, y*1e50, color=:orange, seriestype=:sticks, legend=plotLegend)

    #Hires vibrational scattering
    x = [o2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires]# o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires];
    y = [o2.effCoeff.σ_VibRaman_coeff_0to1_hires]# o2.effCoeff.σ_VibRaman_coeff_1to0_hires];
    @show(x)
    @show(y)
    plot!(x, y*1e50, color=:orange, seriestype=:sticks, legend=plotLegend)

    # Rovibrational scattering (P)
    x = [o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2 o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2];
    y = [o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2 o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2];
    @show(x)
    @show(y)
    #plot!(x, y*1e50, color=:orange, seriestype=:sticks, legend=plotLegend)

    # Rovibrational scattering (R)
    x = [o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2 o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2];
    y = [o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2 o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2];
    @show(x)
    @show(y)
    plot!(x, y*1e50, color=:orange, seriestype=:sticks, legend=plotLegend)

elseif plot_cross_sections==1 
    # N2
    # Rayleigh scattering, pure vibrational Raman scattering
    x = [ν̃+n2.effCoeff.Δν̃_VibRaman_coeff_0to1 ν̃+n2.effCoeff.Δν̃_VibRaman_coeff_1to0 ν̃];
    y = [n2.effCoeff.σ_VibRaman_coeff_0to1*(ν̃+n2.effCoeff.Δν̃_VibRaman_coeff_0to1)^4 n2.effCoeff.σ_VibRaman_coeff_1to0*(ν̃+n2.effCoeff.Δν̃_VibRaman_coeff_1to0)^4 n2.effCoeff.σ_Rayl_coeff*(ν̃)^4];
    @show(ν̃)
    @show(x)
    @show(y)
    #plot((1e7) ./x, y*1e20, color=:blue, seriestype=:sticks, legend=plotLegend)

    # Rayleigh hires scattering
    x = [n2.effCoeff.Δν̃_Rayl_coeff_hires];
    y = [n2.effCoeff.σ_Rayl_coeff_hires];
    @show(ν̃)
    @show(x)
    @show(y)
    plot((1e7) ./x, y*1e50, color=:blue, seriestype=:sticks, legend=plotLegend)

    # Rotational Raman scattering
    x = [ν̃.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2' ν̃.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2'];
    y = [n2.effCoeff.σ_RoRaman_coeff_JtoJp2'.*(ν̃.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2').^4 n2.effCoeff.σ_RoRaman_coeff_JtoJm2'.*(ν̃.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2').^4];
    @show(x)
    @show(y)
    plot!((1e7)*1 ./x, y*1e20, color=:blue, seriestype=:sticks, legend=plotLegend)

    # Rovibrational scattering (P)
    x = [ν̃.+n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2' ν̃.+n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2'];
    y = [n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2'.*(ν̃.+n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2').^4 n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2'.*(ν̃.+n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2').^4];
    @show(x)
    @show(y)
    #plot!((1e7)*1 ./x, y*1e20, color=:blue, seriestype=:sticks, legend=plotLegend)

    # Rovibrational scattering (R)
    x = [ν̃.+n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2' ν̃.+n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2'];
    y = [n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2'.*(ν̃.+n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2').^4 n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2'.*(ν̃.+n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2').^4];
    @show(x)
    @show(y)
    #plot!((1e7)*1 ./x, y*1e20, color=:blue, seriestype=:sticks, legend=plotLegend)


    # O2
    # Rayleigh scattering, pure vibrational Raman scattering
    x = [ν̃+o2.effCoeff.Δν̃_VibRaman_coeff_0to1 ν̃+o2.effCoeff.Δν̃_VibRaman_coeff_1to0 ν̃];
    y = [o2.effCoeff.σ_VibRaman_coeff_0to1*(ν̃+o2.effCoeff.Δν̃_VibRaman_coeff_0to1)^4 o2.effCoeff.σ_VibRaman_coeff_1to0*(ν̃+o2.effCoeff.Δν̃_VibRaman_coeff_1to0)^4 o2.effCoeff.σ_Rayl_coeff*(ν̃)^4];
    #@show(ν̃)
    #@show(x)
    #@show(y)
    #plot!((1e7)*1 ./x, y*1e20, color=:orange, seriestype=:sticks, legend=plotLegend)

    # Rayleigh hires scattering
    x = [o2.effCoeff.Δν̃_Rayl_coeff_hires];
    y = [o2.effCoeff.σ_Rayl_coeff_hires];
    #@show(ν̃)
    #@show(x)
    #@show(y)
    #plot((1e7) ./x, y*1e50, color=:blue, seriestype=:sticks, legend=plotLegend)

    # Rotational Raman scattering
    x = [ν̃.+o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2' ν̃.+o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2'];
    y = [o2.effCoeff.σ_RoRaman_coeff_JtoJp2'.*(ν̃.+o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2').^4 o2.effCoeff.σ_RoRaman_coeff_JtoJm2'.*(ν̃.+o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2').^4];
    #@show(x)
    #@show(y)
    #plot!((1e7)*1 ./x, y*1e20, color=:orange, seriestype=:sticks, legend=plotLegend)

    #Hires vibrational scattering
    x = ν̃.+o2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires# ν̃+o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires];
    y = o2.effCoeff.σ_VibRaman_coeff_0to1_hires.*(ν̃.+o2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires).^4# o2.effCoeff.σ_VibRaman_coeff_1to0_hires];
    @show(x)
    @show(y)
    plot(1e7 ./x, y*1e20, color=:orange, seriestype=:sticks, legend=plotLegend)

    # Rovibrational scattering (P)
    x = [ν̃.+o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2' ν̃.+o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2'];
    y = [o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2'.*(ν̃.+o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2').^4 o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2'.*(ν̃.+o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2').^4];
    @show(x)
    @show(y)
    #plot!((1e7)*1 ./x, y*1e20, color=:orange, seriestype=:sticks, legend=plotLegend)

    # Rovibrational scattering (R)
    x = [ν̃.+o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2' ν̃.+o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2'];
    y = [o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2'.*(ν̃.+o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2').^4 o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2'.*(ν̃.+o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2').^4];
    @show(x)
    @show(y)
    plot!((1e7) ./x, y*1e20, color=:orange, seriestype=:sticks, legend=plotLegend)
end
#foreground_color_legend = nothing
#background_color_legend = nothing