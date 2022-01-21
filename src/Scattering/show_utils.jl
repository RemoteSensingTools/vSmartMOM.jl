#=
 
This file specifies how to pretty-print Scattering module types
 
=#

function Base.show(io::IO, ::MIME"text/plain", x::Aerosol)

    println(io, "\n------------------")
    println(io, "Aerosol Struct")
    println(io, "------------------")
    println(io, "Size distribution: $(x.size_distribution)")
    println(io, "Real refractive Index nᵣ: $(x.nᵣ)")
    println(io, "Imaginary refractive Index nᵣ: $(x.nᵢ)")
    start,stop = quantile(x.size_distribution,[0.0025,0.9975])
    x_axis = collect(range(start,stop,20))
    y = pdf.(x.size_distribution, x_axis) 
    print(io, UnicodePlots.lineplot(x_axis, y,title="Aerosol Size Distribution", xlabel="Radius (μm)" , ylabel="Frequency"))
end

# Overload the show method for ScatteringMatrix
function Base.show(io::IO, x::ScatteringMatrix)
    
    println(io, "  f₁₁| $(length(x.f₁₁))-element vector")
    println(io, "  f₁₂| $(length(x.f₁₂))-element vector")
    println(io, "  f₂₂| $(length(x.f₂₂))-element vector")
    println(io, "  f₃₃| $(length(x.f₃₃))-element vector")
    println(io, "  f₃₄| $(length(x.f₃₄))-element vector")
    print(io, "  f₄₄| $(length(x.f₄₄))-element vector")
end