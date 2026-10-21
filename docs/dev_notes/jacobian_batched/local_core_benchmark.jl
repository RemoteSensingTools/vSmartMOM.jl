# Controlled state-size experiment from supplied core optics: 20 layers,
# three local directions (τ, ϖ, β₂), and 8/20/60 independent input columns.
# All wavelengths batched; full atmospheric adding-doubling and both endpoints.
# Black lower boundary, m=0:2. No optical assembly or upstream physics timed.
# This is an operator experiment, not a microphysical aerosol retrieval scene.
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.Scattering
using YAML, LinearAlgebra, Statistics, Logging
const C = CoreRT
include("core_column_fixture.jl")
const ns = parse(Int,get(ENV,"AUDIT_NSPEC","16"))
const use_gpu = get(ENV,"AUDIT_BACKEND","cpu") == "cuda"
if use_gpu
    using CUDA
    CUDA.functional() || error("CUDA requested but unavailable")
    CUDA.allowscalar(false)
end
sync_device() = use_gpu ? CUDA.synchronize() : nothing
const AT = use_gpu ? CUDA.CuArray : Array
const arch = use_gpu ? GPU() : CPU()
BLAS.set_num_threads(1)
quiet(f) = redirect_stdout(() -> with_logger(f,NullLogger()),devnull)
function timing(f)
    quiet(f); sync_device()
    times = map(1:3) do _
        GC.gc(); sync_device()
        @elapsed begin quiet(f); sync_device() end
    end
    (;median=median(times),samples=times)
end
function main()
    cfg=YAML.load_file(joinpath(pkgdir(vSmartMOM),"test/test_parameters/JacobianTestFast.yaml"))
    delete!(cfg,"absorption"); delete!(cfg,"scattering")
    cfg["radiative_transfer"]["architecture"] = use_gpu ? "GPU()" : "CPU()"
    cfg["radiative_transfer"]["nstreams"] = parse(Int,get(ENV,"AUDIT_STREAMS","3"))
    cfg["radiative_transfer"]["polarization_type"] = get(ENV,"AUDIT_POL","IQU") == "I" ? "Stokes_I()" : "Stokes_IQU()"
    cfg["geometry"]["vaz"] = [0.0,parse(Float64,get(ENV,"AUDIT_AZIMUTH","0"))]
    par=read_parameters(cfg)
    par.spec_bands[1]=collect(range(12987.,13000.;length=ns))
    model=quiet(()->model_from_parameters(par;external_solar=false))
    pol=C.polarization_type(model); q=model.quad_points
    n=length(q.qp_μN); nz=20; FT=Float64
    rs=vSmartMOM.InelasticScattering.noRS{FT}()
    rs.F₀=zeros(FT,pol.n,ns); rs.F₀[1,:].=1
    μ=collect(q.phase_qp_μ)
    g=Scattering.get_greek_rayleigh(0.03)
    zg=Scattering.GreekCoefs(map(f->zeros(FT,length(getfield(g,f))),(:α,:β,:γ,:δ,:ϵ,:ζ))...)
    zg.β[3]=1 # β₀ stays fixed: normalized phase-shape perturbation.
    z=[Scattering.compute_Z_moments(pol,μ,g,m) for m in 0:2]
    dz=[Scattering.compute_Z_moments(pol,μ,zg,m) for m in 0:2]
    ident=Diagonal(AT(ones(FT,n)))
    af=C.make_added_layer(rs,FT,AT,(n,n),ns)
    cf=C.make_composite_layer(rs,FT,AT,(n,n),ns)
    for np in parse.(Int,split(get(ENV,"AUDIT_NPARAMS","8,20,60"),","))
        al,dal=C.make_added_layer(LinMode(),rs,FT,AT,np,(n,n),ns)
        cl,dcl=C.make_composite_layer(LinMode(),rs,FT,AT,np,(n,n),ns)
        local_workspace=C.make_local_jacobian_workspace(rs,FT,AT,3,(n,n),ns,q,pol)
        inputs(delta=zeros(FT,np))=core_optical_inputs(FT,AT,n,ns,nz,np,z,dz,delta)
        data=inputs()
        R=zeros(FT,length(model.obs_geom.vza),pol.n,ns); T=zero(R)
        dR=zeros(FT,size(R)...,np); dT=zero(dR)
        context=(;rs,pol,q,ident,arch,model,af,cf,al,dal,cl,dcl,local_workspace,nz,ns)
        solve(lin,input=data;report=false,factored=false)=
            solve_core_column!(R,T,dR,dT,input,context;linearized=lin,report,factored)
        ft=timing(()->solve(false)); lt=timing(()->solve(true))
        quiet(()->solve(true)); refR=copy(R); refT=copy(T); jac=copy(dR); jacT=copy(dT)
        factored_timing=timing(()->solve(true;factored=true))
        quiet(()->solve(true;factored=true))
        @assert isapprox(dR,jac;rtol=1e-9,atol=1e-11)
        @assert isapprox(dT,jacT;rtol=1e-9,atol=1e-11)
        println("CORE_LOCAL nSpec=$ns layers=$nz pol=$(pol.n) nparams=$np local_directions=3 combined=$factored_timing ratio=$(factored_timing.median/ft.median) speedup=$(lt.median/factored_timing.median)")
        quiet(()->solve(false))
        @assert isapprox(R,refR;rtol=1e-10,atol=1e-12)
        @assert isapprox(T,refT;rtol=1e-10,atol=1e-12)
        println("CORE_MEASURE nSpec=$ns pol=$(pol.n) nparams=$np forward=$ft combined=$lt ratio=$(lt.median/ft.median)"); flush(stdout)
        if get(ENV,"AUDIT_CORE_TIMERS","false") == "true"
            println("CORE_PROFILE nparams=$np forward")
            with_logger(()->solve(false;report=true),NullLogger())
            println("CORE_PROFILE nparams=$np combined")
            with_logger(()->solve(true;report=true),NullLogger())
        end
        if ns<=64 && get(ENV,"AUDIT_FD","true")=="true"
            for p in 1:np
                δ=zeros(np); δ[p]=1e-5
                plus,minus=inputs(δ),inputs(-δ)
                quiet(()->solve(false,plus)); rp=copy(R); tp=copy(T)
                quiet(()->solve(false,minus)); rm=copy(R); tm=copy(T)
                fd=(rp-rm)/(2δ[p])
                @assert isapprox(fd,jac[:,:,:,p];rtol=2e-5,atol=2e-8)
                @assert isapprox((tp-tm)/(2δ[p]),jacT[:,:,:,p];rtol=2e-5,atol=2e-8)
            end
            println("CORE_FD_PASS nparams=$np all_columns=true endpoints=R,T azimuth=$(model.obs_geom.vaz)")
        end
    end
end
main()
