using Plots

γRayl_O2 = 0.02885
γRayl_N2 = 0.01058
γVib_O2 = 0.16951
γVib_N2 = 0.08497
γRot = 0.75

θ = 0:1:360

F11_Rayl_O2 = (3/(8π*(2γRayl_O2+1)))*((cosd.(θ)).^2 .+ γRayl_O2*(sind.(θ)).^2)
F11_Rayl_N2 = (3/(8π*(2γRayl_N2+1)))*((cosd.(θ)).^2 .+ γRayl_N2*(sind.(θ)).^2)
F11_Vib_O2  = (3/(8π*(2γVib_O2+1)))*((cosd.(θ)).^2 .+ γVib_O2*(sind.(θ)).^2)
F11_Vib_N2  = (3/(8π*(2γVib_N2+1)))*((cosd.(θ)).^2 .+ γVib_N2*(sind.(θ)).^2)
F11_Rot     = (3/(8π*(2γRot+1)))*((cosd.(θ)).^2 .+ γRot*(sind.(θ)).^2)

F22_Rayl_O2 = (3/(8π*(2γRayl_O2+1)))*ones(length(θ))
F22_Rayl_N2 = (3/(8π*(2γRayl_N2+1)))*ones(length(θ))
F22_Vib_O2  = (3/(8π*(2γVib_O2+1)))*ones(length(θ))
F22_Vib_N2  = (3/(8π*(2γVib_N2+1)))*ones(length(θ))
F22_Rot     = (3/(8π*(2γRot+1)))*ones(length(θ))

F12_Rayl_O2 = (3/(8π*(2γRayl_O2+1)))*γRayl_O2*ones(length(θ))
F12_Rayl_N2 = (3/(8π*(2γRayl_N2+1)))*γRayl_N2*ones(length(θ))
F12_Vib_O2  = (3/(8π*(2γVib_O2+1)))*γVib_O2*ones(length(θ))
F12_Vib_N2  = (3/(8π*(2γVib_N2+1)))*γVib_N2*ones(length(θ))
F12_Rot     = (3/(8π*(2γRot+1)))*γRot*ones(length(θ))

F33_Rayl_O2 = (3/(8π*(2γRayl_O2+1)))*(1-γRayl_O2).*cosd.(θ)
F33_Rayl_N2 = (3/(8π*(2γRayl_N2+1)))*(1-γRayl_N2).*cosd.(θ)
F33_Vib_O2  = (3/(8π*(2γVib_O2+1)))*(1-γVib_O2).*cosd.(θ)
F33_Vib_N2  = (3/(8π*(2γVib_N2+1)))*(1-γVib_N2).*cosd.(θ)
F33_Rot     = (3/(8π*(2γRot+1)))*(1-γRot).*cosd.(θ)

F44_Rayl_O2 = (3/(8π*(2γRayl_O2+1)))*(1-3γRayl_O2).*cosd.(θ)
F44_Rayl_N2 = (3/(8π*(2γRayl_N2+1)))*(1-3γRayl_N2).*cosd.(θ)
F44_Vib_O2  = (3/(8π*(2γVib_O2+1)))*(1-3γVib_O2).*cosd.(θ)
F44_Vib_N2  = (3/(8π*(2γVib_N2+1)))*(1-3γVib_N2).*cosd.(θ)
F44_Rot     = (3/(8π*(2γRot+1)))*(1-3γRot).*cosd.(θ)


l = @layout [a1 a2]

f1 = plot(θ, F11_Rayl_O2, color=:red, label="O₂ Rayl", linewidth=3, xlabel="θ [ᵒ]")
f1 = plot!(θ, F11_Rayl_N2, color=:blue, label="N₂ Rayl", linewidth=3)
f1 = plot!(θ, F11_Vib_O2, color=:red, linestyle=:dot, label="O₂ VRS", linewidth=3)
f1 = plot!(θ, F11_Vib_N2, color=:blue, linestyle=:dot, label="N₂ VRS", linewidth=3)
f1 = plot!(θ, F11_Rot, color=:black, label="RRS/RVRS", linewidth=3)

f2 = plot(θ*π/180, F11_Rayl_O2, proj=:polar, color=:red, label="", linewidth=3)
f2 = plot!(θ*π/180, F11_Rayl_N2, proj=:polar, color=:blue, label="", linewidth=3)
f2 = plot!(θ*π/180, F11_Vib_O2, proj=:polar, color=:red, label="", linestyle=:dot, linewidth=3)
f2 = plot!(θ*π/180, F11_Vib_N2, proj=:polar, color=:blue, label="", linestyle=:dot, linewidth=3)
f2 = plot!(θ*π/180, F11_Rot, proj=:polar, color=:black, label="", linewidth=3)
plot(f1, f2, layout = l, title = ["F₁₁" "F₁₁"], foreground_color_legend = nothing, background_color_legend = nothing, titlefont = font(10))
savefig("RamanPhaseF11.png")

l = @layout [a1 a2]

f1 = plot(θ, F22_Rayl_O2, color=:red, label="O₂ Rayl", linewidth=3, xlabel="θ [ᵒ]")
f1 = plot!(θ, F22_Rayl_N2, color=:blue, label="N₂ Rayl", linewidth=3)
f1 = plot!(θ, F22_Vib_O2, color=:red, linestyle=:dot, label="O₂ VRS", linewidth=3)
f1 = plot!(θ, F22_Vib_N2, color=:blue, linestyle=:dot, label="N₂ VRS", linewidth=3)
f1 = plot!(θ, F22_Rot, color=:black, label="RRS/RVRS", linewidth=3)

f2 = plot(θ*π/180, F22_Rayl_O2, proj=:polar, color=:red, label="", linewidth=3)
f2 = plot!(θ*π/180, F22_Rayl_N2, proj=:polar, color=:blue, label="", linewidth=3)
f2 = plot!(θ*π/180, F22_Vib_O2, proj=:polar, color=:red, label="", linestyle=:dot, linewidth=3)
f2 = plot!(θ*π/180, F22_Vib_N2, proj=:polar, color=:blue, label="", linestyle=:dot, linewidth=3)
f2 = plot!(θ*π/180, F22_Rot, proj=:polar, color=:black, label="", linewidth=3)
plot(f1, f2, layout = l, title = ["F₂₂" "F₂₂"], foreground_color_legend = nothing, background_color_legend = nothing, titlefont = font(10))
savefig("RamanPhaseF22.png")

l = @layout [a1 a2]

f1 = plot(θ, F12_Rayl_O2, color=:red, label="O₂ Rayl", linewidth=3, xlabel="θ [ᵒ]")
f1 = plot!(θ, F12_Rayl_N2, color=:blue, label="N₂ Rayl", linewidth=3)
f1 = plot!(θ, F12_Vib_O2, color=:red, linestyle=:dot, label="O₂ VRS", linewidth=3)
f1 = plot!(θ, F12_Vib_N2, color=:blue, linestyle=:dot, label="N₂ VRS", linewidth=3)
f1 = plot!(θ, F12_Rot, color=:black, label="RRS/RVRS", linewidth=3)

f2 = plot(θ*π/180, F12_Rayl_O2, proj=:polar, color=:red, label="", linewidth=3)
f2 = plot!(θ*π/180, F12_Rayl_N2, proj=:polar, color=:blue, label="", linewidth=3)
f2 = plot!(θ*π/180, F12_Vib_O2, proj=:polar, color=:red, linestyle=:dot, label="", linewidth=3)
f2 = plot!(θ*π/180, F12_Vib_N2, proj=:polar, color=:blue, linestyle=:dot, label="", linewidth=3)
f2 = plot!(θ*π/180, F12_Rot, proj=:polar, color=:black, label="", linewidth=3)
plot(f1, f2, layout = l, title = ["F₁₂" "F₁₂"], foreground_color_legend = nothing, background_color_legend = nothing, titlefont = font(10))

savefig("RamanPhaseF12.png")

l = @layout [a1 a2]

f1 = plot(θ, F33_Rayl_O2, color=:red, label="O₂ Rayl", linewidth=3, xlabel="θ [ᵒ]")
f1 = plot!(θ, F33_Rayl_N2, color=:blue, label="N₂ Rayl", linewidth=3)
f1 = plot!(θ, F33_Vib_O2, color=:red, linestyle=:dot, label="O₂ VRS", linewidth=3)
f1 = plot!(θ, F33_Vib_N2, color=:blue, linestyle=:dot, label="N₂ VRS", linewidth=3)
f1 = plot!(θ, F33_Rot, color=:black, label="RRS/RVRS", linewidth=3)

f2 = plot(θ*π/180, F33_Rayl_O2, proj=:polar, color=:red, label="", linewidth=3)
f2 = plot!(θ*π/180, F33_Rayl_N2, proj=:polar, color=:blue, label="", linewidth=3)
f2 = plot!(θ*π/180, F33_Vib_O2, proj=:polar, color=:red, linestyle=:dot, label="", linewidth=3)
f2 = plot!(θ*π/180, F33_Vib_N2, proj=:polar, color=:blue, linestyle=:dot, label="", linewidth=3)
f2 = plot!(θ*π/180, F33_Rot, proj=:polar, color=:black, label="", linewidth=3)
f2 = plot!(θ*π/180, -F33_Rayl_O2, proj=:polar, color=:orange, label="", linewidth=3)
f2 = plot!(θ*π/180, -F33_Rayl_N2, proj=:polar, color=:green, label="", linewidth=3)
f2 = plot!(θ*π/180, -F33_Vib_O2, proj=:polar, color=:orange, label="", linestyle=:dot, linewidth=3)
f2 = plot!(θ*π/180, -F33_Vib_N2, proj=:polar, color=:green, label="", linestyle=:dot, linewidth=3)
f2 = plot!(θ*π/180, -F33_Rot, proj=:polar, color=:brown, label="", linewidth=3)

plot(f1, f2, layout = l, title = ["F₃₃" "F₃₃"], foreground_color_legend = nothing, background_color_legend = nothing, titlefont = font(10))

savefig("RamanPhaseF33.png")

l = @layout [a1 a2]

f1 = plot(θ, F44_Rayl_O2, color=:red, label="O₂ Rayl", linewidth=3)
f1 = plot!(θ, F44_Rayl_N2, color=:blue, label="N₂ Rayl", linewidth=3)
f1 = plot!(θ, F44_Vib_O2, color=:red, linestyle=:dot, label="O₂ VRS", linewidth=3)
f1 = plot!(θ, F44_Vib_N2, color=:blue, linestyle=:dot, label="N₂ VRS", linewidth=3)
f1 = plot!(θ, F44_Rot, color=:black, label="RRS/RVRS", linewidth=3)

f2 = plot(θ*π/180, F44_Rayl_O2, proj=:polar, color=:red, label="", linewidth=3)
f2 = plot!(θ*π/180, F44_Rayl_N2, proj=:polar, color=:blue, label="", linewidth=3)
f2 = plot!(θ*π/180, F44_Vib_O2, proj=:polar, color=:red, linestyle=:dot, label="", linewidth=3)
f2 = plot!(θ*π/180, F44_Vib_N2, proj=:polar, color=:blue, linestyle=:dot, label="", linewidth=3)
f2 = plot!(θ*π/180, F44_Rot, proj=:polar, color=:black, label="", linewidth=3)
f2 = plot!(θ*π/180, -F44_Rayl_O2, proj=:polar, color=:orange, label="", linewidth=3)
f2 = plot!(θ*π/180, -F44_Rayl_N2, proj=:polar, color=:green, label="", linewidth=3)
f2 = plot!(θ*π/180, -F44_Vib_O2, proj=:polar, color=:orange, linestyle=:dot, label="", linewidth=3)
f2 = plot!(θ*π/180, -F44_Vib_N2, proj=:polar, color=:green, linestyle=:dot, label="", linewidth=3)
f2 = plot!(θ*π/180, -F44_Rot, proj=:polar, color=:brown, label="", linewidth=3)

plot(f1, f2, layout = l, title = ["F₄₄" "F₄₄"], foreground_color_legend = nothing, background_color_legend = nothing, titlefont = font(10))

savefig("RamanPhaseF44.png")
#plot(f11, f12, f33, f44, layout = l, legend = false, title = ["F₁₁" "F₁₂" "F₃₃" "F₄₄"], titlefont = font(10))

