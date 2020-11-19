module RTM
using LinearAlgebra
using ..PhaseFunction

using FastGaussQuadrature 

include("types.jl") 
include("RTM_main.jl")
include("atmo_prof.jl")
include("rt_elemental.jl")
include("rt_doubling.jl")
include("rt_interaction.jl")
include("rt_tools.jl")
#include("")
#include("RTM_main.jl")
#include("RTM_main.jl")

export rt_set_streams

end