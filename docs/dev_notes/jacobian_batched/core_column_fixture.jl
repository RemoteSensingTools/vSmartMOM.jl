# Shared controlled core-optics fixture for the local-basis and layer-response
# experiments. Every wavelength repeats the same optics; the benchmark isolates
# batch size and retrieval-column count, not spectroscopic variability.

"""
    core_optical_inputs(FT, AT, n, ns, nz, np, z, dz, delta=zeros(FT,np))

Build τ, ϖ and normalized β₂ phase perturbations in every layer. A cyclic
column map makes retrieval parameters affect multiple layers when np < 3nz.
Return the same perturbation both as complete physical tangents and as three
local directions with scalar chain coefficients. `delta` perturbs the forward
inputs independently for finite-difference checks; β₀ remains fixed.
"""
function core_optical_inputs(FT,AT,n,ns,nz,np,z,dz,delta=zeros(FT,np))
    values=[]; tangents=[]; local_tangents=[]
    for m in 0:2
        vm=[]; dm=[]; lm=[]
        for k in 1:nz
            dt=zeros(FT,ns,np); dw=zero(dt)
            dp=zeros(FT,n,n,1,np); dn=zero(dp)
            columns = [mod(3(k-1)+j-1,np)+1 for j in 1:3]
            dt[:,columns[1]].=1
            dw[:,columns[2]].=1
            dp[:,:,1,columns[3]].=dz[m+1][1]
            dn[:,:,1,columns[3]].=dz[m+1][2]
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
            bt=zeros(FT,ns,3); bw=zero(bt); bt[:,1].=1; bw[:,2].=1
            bp=zeros(FT,n,n,1,3); bn=zero(bp)
            bp[:,:,:,3].=reshape(dz[m+1][1],n,n,1)
            bn[:,:,:,3].=reshape(dz[m+1][2],n,n,1)
            coeff=zeros(FT,ns,3,np)
            for j in 1:3
                coeff[:,j,columns[j]].=1
            end
            basis=C.CoreScatteringOpticalPropertiesLin(AT(bt),AT(bw),AT(bp),AT(bn))
            push!(lm,C.LocalOpticalJacobian(AT(dt),AT(dw),basis,AT(coeff)))
        end
        push!(values,vm); push!(tangents,dm); push!(local_tangents,lm)
    end
    interfaces,sums,dsums=C.extractEffectiveProps(values[1],tangents[1])
    (;values,tangents,local_tangents,interfaces,sums=[sums[:,k] for k in 1:nz],
      dsums=[dsums[:,:,k] for k in 1:nz],
      maxima=[maximum(v.τ.*v.ϖ) for v in values[1]])
end

"""
    solve_core_column!(R, T, dR, dT, data, context; linearized=false,
                       report=false, factored=false)

Run m=0:2 over the supplied core column and accumulate both endpoint radiances.
The lower boundary is black; the forward and tangent paths share quadrature,
solar source and numerical settings. `factored=true` changes only the local
optical handoff to doubling. Preparation and workspace allocation occur outside
this function so both benchmark drivers time the same complete core solve.
"""
function solve_core_column!(R,T,dR,dT,data,context;linearized=false,report=false,factored=false)
    (;rs,pol,q,ident,arch,model,af,cf,al,dal,cl,dcl,local_workspace,nz,ns)=context
    lin=linearized
    fill!(R,0); fill!(T,0); lin && (fill!(dR,0);fill!(dT,0))
    for m in 0:2
        for k in 1:nz
            if lin
                C.rt_kernel!(rs,pol,true,al,dal,cl,dcl,
                    data.values[m+1][k],(factored ? data.local_tangents : data.tangents)[m+1][k],data.interfaces[k],
                    data.sums[k],data.dsums[k],m,q,ident,arch,q.qp_μN,k;
                    local_workspace=factored ? local_workspace : nothing)
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
