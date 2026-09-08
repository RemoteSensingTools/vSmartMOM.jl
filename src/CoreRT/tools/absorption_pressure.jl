"""
    _absorption_pressure_derivative(model, grid, p, T, broadener)

Cross-section derivative with respect to layer pressure [hPa], at fixed
temperature and composition. This upstream boundary returns ordinary arrays;
no dual numbers enter the radiative-transfer operator.
"""
function _absorption_pressure_derivative(model, grid, p, T, broadener)
    throw(ArgumentError("pressure derivatives are not implemented for $(typeof(model))"))
end

# Legacy BSpline LUTs accept dual pressure coordinates. Differentiate their
# actual interpolant, including its spectral out-of-range behavior.
function _absorption_pressure_derivative(
    model::Absorption.InterpolationModel, grid, p, T, broadener)
    pressure = ForwardDiff.Dual{Nothing}(p, one(p))
    σ = _layer_absorption_cross_section(model, grid, pressure, T, broadener)
    return ForwardDiff.partials.(σ, 1)
end

# ABSCO uses a linear pressure blend of spectra evaluated at each pressure
# node's OWN temperature grid. Evaluating the two endpoints through the public
# API preserves those temperature, broadener, and spectral interpolation rules.
# At an interior knot choose the right interval; at/outside either endpoint
# choose the clamped (zero) derivative. A centered difference at a kink need
# not agree with this explicitly one-sided convention.
function _absorption_pressure_derivative(
    model::Union{AtmosphericAbsorption.AbscoLUT{FT},
                 AtmosphericAbsorption.InterpolationModel{FT}},
    grid, p, T, broadener) where FT
    pressure = FT(p)
    nodes = model.p
    if pressure <= first(nodes) || pressure >= last(nodes)
        return zeros(FT, length(grid))
    end
    lo = searchsortedlast(nodes, pressure)
    hi = lo + 1
    σlo = _layer_absorption_cross_section(model, grid, nodes[lo], T, broadener)
    σhi = _layer_absorption_cross_section(model, grid, nodes[hi], T, broadener)
    return Array((σhi .- σlo) ./ (nodes[hi] - nodes[lo]))
end

# AA prepares pressure-linear collisional parameters at fixed T and broadener.
# Seed a dimensionless pressure scale: dΓ/dscale=Γ, dΔ/dscale=Δ, etc.
# Differentiating the actual profile/CPF implementation keeps its approximation
# consistent with the forward spectrum. Divide by p to obtain dσ/dp [hPa⁻¹].

# The optional CPU reference CPF calls complex erfcx, which SpecialFunctions
# does not accept on Complex{Dual}. Supply its exact scalar chain rule through
# our own strategy type; do not add methods to foreign numeric types.
struct _PressureErfcxCPF <: AtmosphericAbsorption.LineShapes.AbstractCPF end
_pressure_cpf(cpf) = cpf
_pressure_cpf(::AtmosphericAbsorption.ErfcxCPF) = _PressureErfcxCPF()

@inline function AtmosphericAbsorption.LineShapes.w(::_PressureErfcxCPF,
    z::Complex{ForwardDiff.Dual{Tag,FT,N}}) where {Tag,FT,N}
    z₀ = complex(ForwardDiff.value(real(z)), ForwardDiff.value(imag(z)))
    w₀ = AtmosphericAbsorption.LineShapes.w(AtmosphericAbsorption.ErfcxCPF(), z₀)
    # w(z)=erfcx(-iz), so w′(z)=-2z*w(z)+2i/√π.
    slope = -2z₀ * w₀ + complex(zero(FT), FT(2) / sqrt(FT(π)))
    dz = ntuple(k -> complex(ForwardDiff.partials(real(z), k),
                             ForwardDiff.partials(imag(z), k)), Val(N))
    real_dot = ForwardDiff.Partials(ntuple(k -> real(slope * dz[k]), Val(N)))
    imag_dot = ForwardDiff.Partials(ntuple(k -> imag(slope * dz[k]), Val(N)))
    return complex(ForwardDiff.Dual{Tag}(real(w₀), real_dot),
                   ForwardDiff.Dual{Tag}(imag(w₀), imag_dot))
end

@inline function _line_pressure_derivative(profile, cpf, ν, ν₀, γd,
                                          Γ0, Γ2, Δ0, Δ2, νVC, η, Y)
    FT = typeof(ν)
    fixed(x) = ForwardDiff.Dual{Nothing}(x, zero(FT))
    scaled(x) = ForwardDiff.Dual{Nothing}(x, x)
    parameters = (γd=fixed(γd), Γ0=scaled(Γ0), Γ2=scaled(Γ2),
                  Δ0=scaled(Δ0), Δ2=scaled(Δ2), νVC=scaled(νVC),
                  η=fixed(η), Y=scaled(Y))
    value = AtmosphericAbsorption.LineShapes.evaluate(
        profile, _pressure_cpf(cpf), fixed(ν), fixed(ν₀), parameters)
    return ForwardDiff.partials(value, 1)
end

@kernel function _line_pressure_kernel!(σp, @Const(grid), @Const(ν₀), @Const(γd),
    @Const(Γ0), @Const(Γ2), @Const(Δ0), @Const(Δ2), @Const(νVC), @Const(η),
    @Const(Y), @Const(S), @Const(istart), @Const(istop), n, profile, cpf, p)
    i = @index(Global, Linear)
    value = zero(eltype(σp))
    @inbounds for j in 1:n
        if istart[j] <= i <= istop[j]
            value += S[j] * _line_pressure_derivative(
                profile, cpf, grid[i], ν₀[j], γd[j], Γ0[j], Γ2[j],
                Δ0[j], Δ2[j], νVC[j], η[j], Y[j])
        end
    end
    @inbounds σp[i] = value / p
end

function _absorption_pressure_derivative(
    model::AtmosphericAbsorption.LineByLineModel{FT}, grid, p, T, broadener) where FT
    pressure = FT(p)
    pressure > zero(FT) || throw(ArgumentError("line pressure must be positive"))
    isempty(grid) && return FT[]
    vmr = broadener === nothing ? model.vmr : broadener
    # Keep the active-line windows from the forward preparation. Hard wing
    # cutoffs are discontinuous; their membership has no derivative at a jump.
    prepared = AtmosphericAbsorption.Crosssections.prepare(model, grid, p, T; vmr)
    AAArch = AtmosphericAbsorption.Architectures
    backend = AAArch.devi(model.architecture)
    to = AAArch.array_type(model.architecture)
    σp = to(zeros(FT, length(grid)))
    _line_pressure_kernel!(backend)(σp, to(collect(FT, grid)),
        prepared.ν0, prepared.γd, prepared.Γ0, prepared.Γ2, prepared.Δ0,
        prepared.Δ2, prepared.νVC, prepared.η, prepared.Y, prepared.S,
        prepared.istart, prepared.istop, prepared.n, model.profile, model.cpf,
        pressure; ndrange=length(grid))
    KernelAbstractions.synchronize(backend)
    return Array(σp)
end

_accumulate_absorption_pressure!(::Nothing, model, grid, vmr, profile, broadener) = nothing

"Accumulate only the cross-section part of dτ/dp_surf; column scaling is added by the constructor."
function _accumulate_absorption_pressure!(τp, model, grid, vmr, profile, broadener)
    iz = length(profile.p_full)
    abundance = vmr isa AbstractArray ? vmr[iz] : vmr
    iszero(abundance) && return nothing
    broadener_layer = broadener isa AbstractArray ? broadener[iz] : broadener
    broadener_layer = broadener_layer === nothing ?
        _default_broadener(model, profile, iz) : broadener_layer
    σp = _absorption_pressure_derivative(
        model, grid, profile.p_full[iz], profile.T[iz], broadener_layer)
    # Fixed final grid: only p_half[end] varies, hence dp_full[end]/dp_surf=1/2.
    column = profile.vcd_dry[iz] * abundance / 2
    view(τp, :, iz) .+= σp .* column
    return nothing
end
