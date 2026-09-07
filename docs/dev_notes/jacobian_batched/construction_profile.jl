# Run from test/ with the same AUDIT_* scene settings as the RT benchmark.
# Load its fixture/helpers without launching its RT timing loop.
include(ex -> ex == :(main()) ? nothing : ex,
        joinpath(@__DIR__, "local_basis_benchmark.jl"))

println("AtmosphericAbsorption source: ", pathof(CoreRT.AtmosphericAbsorption))
for mode in (:forward, :linearized)
    construct = () -> begin
        p = fixture_parameters(1)
        mode === :forward ? model_from_parameters(p; external_solar) :
            model_from_parameters(LinMode(),p; external_solar)
    end
    quiet(construct); sync_backend()
    for sample in 1:3
        GC.gc(); sync_backend()
        CoreRT.reset_timer!()
        t = @timed begin quiet(construct); sync_backend() end
        println("CONSTRUCTION mode=$mode sample=$sample seconds=$(t.time) host_bytes=$(t.bytes) gc_seconds=$(t.gctime)")
        CoreRT.print_timer()
        println(); flush(stdout)
    end
    Profile.clear()
    @profile begin quiet(construct); sync_backend() end
    println("HOST_PROFILE mode=$mode")
    Profile.print(stdout; format=:flat, C=false, sortedby=:count, mincount=10)
    flush(stdout)
end
