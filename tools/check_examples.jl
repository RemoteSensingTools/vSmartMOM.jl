using vSmartMOM, Test
text = read(joinpath(pkgdir(vSmartMOM), "README.md"), String)
@testset "Executable README RT examples" begin
    for heading in ("#### Forward run (minimal)", "#### Linearized run (analytic Jacobians)")
        section = split(text, heading; limit=2)[2]
        code = match(r"```julia\n(.*?)```"s, section).captures[1]
        mod = Module(gensym(:ReadmeExample))
        Base.include_string(mod, code, "README.md")
        @test all(isfinite, Base.invokelatest(getfield, mod, :R))
        @test all(isfinite, Base.invokelatest(getfield, mod, :T))
        if occursin("Linearized", heading)
            @test all(isfinite, Base.invokelatest(getfield, mod, :dR))
            @test all(isfinite, Base.invokelatest(getfield, mod, :dT))
            @test size(Base.invokelatest(getfield, mod, :dR), 4) == 1
        end
    end
end
