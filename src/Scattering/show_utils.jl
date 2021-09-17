#=
 
This file specifies how to pretty-print Scattering module types
 
=#

# Overload the show method for ScatteringMatrix
function Base.show(io::IO, x::ScatteringMatrix)
    
    println(io, "  f₁₁| $(length(x.f₁₁))-element vector")
    println(io, "  f₁₂| $(length(x.f₁₂))-element vector")
    println(io, "  f₂₂| $(length(x.f₂₂))-element vector")
    println(io, "  f₃₃| $(length(x.f₃₃))-element vector")
    println(io, "  f₃₄| $(length(x.f₃₄))-element vector")
    print(io, "  f₄₄| $(length(x.f₄₄))-element vector")
end