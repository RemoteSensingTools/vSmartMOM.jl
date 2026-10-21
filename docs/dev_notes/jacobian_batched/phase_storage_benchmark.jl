# Separate angular node evaluation, spectral expansion and layer mixing.
# Run from test/: CUDA_VISIBLE_DEVICES=0 AUDIT_NSPEC=10000 julia --project=. ...
using vSmartMOM, vSmartMOM.CoreRT, Test, YAML, Logging, LinearAlgebra, Statistics, TOML, CUDA
include(joinpath(pkgdir(vSmartMOM),"test/local_jacobian_fixture.jl"))
CUDA.allowscalar(false)
BLAS.set_num_threads(1)
const C = CoreRT
quiet(f) = with_logger(f,NullLogger())
ns = parse(Int,get(ENV,"AUDIT_NSPEC","10000"))
nz = parse(Int,get(ENV,"AUDIT_LAYERS","20"))
external = get(ENV,"AUDIT_EXTERNAL_SOLAR","false") == "true"
model,lin = quiet(()->local_jacobian_fixture(Float64,true,external;
    gpu=true,nspec=ns,nlayers=nz))
AT = C.array_type(model)
rs = vSmartMOM.InelasticScattering.noRS{Float64}()
cache = C.build_local_jacobian_cache(1,model,lin,C.build_m_invariant_cache_lin(1,model,lin))
m = 1
optics, tangent = only(model.aerosol_optics[1]), only(lin.lin_aerosol_optics[1])
ray = C._compute_phase_blocks(model,model.greek_rayleigh[1],m,AT)
aerosol = C._compute_aerosol_phase_blocks_lin(model,optics,tangent,C.get_spec_bands(model)[1],m,AT)
make_basis() = ntuple(k->C.local_phase_basis(ray[k],
    [(aerosol[(1,2,5,6)[k]],aerosol[(3,4,7,8)[k]])]),4)
basis = make_basis()
forward_phases = map(k->[aerosol[k]],(1,2,5,6))
function measure_phase(f)
    f(); CUDA.synchronize()
    times, allocations = Float64[], Int[]
    for _ in 1:5
        GC.gc(); CUDA.synchronize()
        result = CUDA.@timed f()
        push!(times,result.time)
        push!(allocations,result.gpu_bytes)
    end
    Dict("seconds"=>times,"median_seconds"=>median(times),"allocated_device_bytes"=>allocations)
end
stages = Dict(
    "angular_nodes"=>()->[C._compute_phase_blocks_lin(model,g,dg,m,AT)
        for (g,dg) in zip(optics.phase_greek,tangent.phase_lin_greek)],
    "angular_nodes_and_spectral_expansion"=>()->C._compute_aerosol_phase_blocks_lin(
        model,optics,tangent,C.get_spec_bands(model)[1],m,AT),
    "local_basis"=>make_basis,
    "all_layer_mixtures"=>()->[ntuple(k->C.mix_local_forward_phase(
        ray[k],forward_phases[k],layer.mixing_weights),4)
        for layer in cache.layers],
    "complete_phase_assembly"=>()->C.construct_local_optical_jacobians(rs,1,m,model,lin,cache))
records = Dict(name=>measure_phase(f) for (name,f) in stages)
properties,jacobians,_ = C.construct_local_optical_jacobians(rs,1,m,model,lin,cache)
@assert all(j->j.basis === first(jacobians).basis,jacobians)
payload(a) = a === nothing ? 0 : sizeof(eltype(a))*length(a)
result = Dict("spectral_points"=>ns,"layers"=>nz,"external_solar"=>external,
    "aerosol_nodes"=>length(optics.phase_ν),"fourier_order"=>m,
    "shared_phase_basis_bytes"=>sum(payload,basis),
    "all_layer_mixture_bytes"=>sum(p->sum(payload,(p.Z⁺⁺,p.Z⁻⁺,p.Z₀⁺,p.Z₀⁻)),properties),
    "stages"=>records)
println(result)
open(io->TOML.print(io,result),get(ENV,"AUDIT_OUTPUT","/tmp/phase-storage.toml"),"w")
