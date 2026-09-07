# Full atmospheric adding–doubling from supplied core properties.
# Run from test/: AUDIT_BACKEND=cuda AUDIT_NSPEC=10000 AUDIT_POL=IQU julia --project=. ../docs/dev_notes/jacobian_batched/core_properties_benchmark.jl
# Black lower boundary; solar SFI, m=0:2, TOA/BOA diffuse postprocessing.
# No Mie, mixing, or physical-parameter expansion is timed.
# Directions: τ, ϖ, degree-2 Greek β. β₂ is ONE normalized phase-shape
# direction, not the dense dI/dZ tensor. np=15 varies these in each layer.
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.Scattering
using YAML, LinearAlgebra, Statistics, Logging
const C = CoreRT
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
    cfg["radiative_transfer"]["polarization_type"] = get(ENV,"AUDIT_POL","IQU") == "I" ? "Stokes_I()" : "Stokes_IQU()"
    cfg["geometry"]["vaz"] = [0.0,parse(Float64,get(ENV,"AUDIT_AZIMUTH","0"))]
    par=read_parameters(cfg)
    par.spec_bands[1]=collect(range(12987.,13000.;length=ns))
    model=quiet(()->model_from_parameters(par;external_solar=false))
    pol=C.polarization_type(model); q=model.quad_points
    n=length(q.qp_μN); nz=5; FT=Float64
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
    for np in parse.(Int,split(get(ENV,"AUDIT_NPARAMS","1,3,15"),","))
        al,dal=C.make_added_layer(LinMode(),rs,FT,AT,np,(n,n),ns)
        cl,dcl=C.make_composite_layer(LinMode(),rs,FT,AT,np,(n,n),ns)
        function inputs(delta=zeros(np))
            values=[]; tangents=[]
            for m in 0:2
                vm=[]; dm=[]
                for k in 1:nz
                    dt=zeros(FT,ns,np); dw=zero(dt)
                    dp=zeros(FT,n,n,1,np); dn=zero(dp)
                    base=np==15 ? 3(k-1) : 0
                    dt[:,base+1].=1
                    if np>1
                        dw[:,base+2].=1
                        dp[:,:,1,base+3].=dz[m+1][1]
                        dn[:,:,1,base+3].=dz[m+1][2]
                    end
                    τ=fill(0.025+0.012k,ns)+dt*delta
                    ω=fill(0.9,ns)+dw*delta
                    zp=reshape(copy(z[m+1][1]),n,n,1)
                    zn=reshape(copy(z[m+1][2]),n,n,1)
                    for p in 1:np
                        zp .+= delta[p].*dp[:,:,:,p]
                        zn .+= delta[p].*dn[:,:,:,p]
                    end
                    push!(vm,C.CoreScatteringOpticalProperties(AT(τ),AT(ω),AT(zp),AT(zn)))
                    push!(dm,C.CoreScatteringOpticalPropertiesLin(AT(dt),AT(dw),AT(dp),AT(dn)))
                end
                push!(values,vm); push!(tangents,dm)
            end
            interfaces,sums,dsums=C.extractEffectiveProps(values[1],tangents[1])
            (;values,tangents,interfaces,sums=[sums[:,k] for k in 1:nz],
              dsums=[dsums[:,:,k] for k in 1:nz],
              maxima=[maximum(v.τ.*v.ϖ) for v in values[1]])
        end
        data=inputs()
        R=zeros(FT,length(model.obs_geom.vza),pol.n,ns); T=zero(R)
        dR=zeros(FT,size(R)...,np); dT=zero(dR)
        function solve(lin,data=data; report=false)
            fill!(R,0); fill!(T,0); lin && (fill!(dR,0);fill!(dT,0))
            for m in 0:2
                for k in 1:nz
                    if lin
                        C.rt_kernel!(rs,pol,true,al,dal,cl,dcl,
                            data.values[m+1][k],data.tangents[m+1][k],data.interfaces[k],
                            data.sums[k],data.dsums[k],m,q,ident,arch,q.qp_μN,k)
                    else
                        C.rt_kernel!(rs,pol,true,af,cf,data.values[m+1][k],
                            data.interfaces[k],data.sums[k],m,q,ident,arch,q.qp_μN,k;
                            max_τϖ=data.maxima[k])
                    end
                end
                args=(model.obs_geom.vza,q.qp_μ,m,model.obs_geom.vaz,q.μ₀,
                      m==0 ? 0.5/π : 1/π,ns,true)
                if lin
                    C.postprocessing_vza!(rs,q.iμ₀,pol,cl,dcl,args...,R,T,dR,dT)
                else
                    C.postprocessing_vza!(rs,q.iμ₀,pol,cf,args...,nothing,R,nothing,T,nothing,nothing)
                end
            end
            report && C.print_timer()
            C.reset_timer!()
            nothing
        end
        ft=timing(()->solve(false)); lt=timing(()->solve(true))
        quiet(()->solve(true)); refR=copy(R); refT=copy(T); jac=copy(dR); jacT=copy(dT)
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
        if ns<=64
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
