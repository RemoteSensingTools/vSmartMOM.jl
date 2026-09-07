# Private reference switch for numerical and allocation comparisons.
const _Z_JACOBIAN_TABLES_ENABLED = Ref(true)

# The phase operator is linear in the Greek coefficients. Forward values and
# their tangents therefore share the same angular tables and static products.
# Sanghavi, JQSRT 136 (2014), 16–27, (14)–(16):
# Aᵐ(μᵢ,μⱼ) = Σₗ Πᵐₗ(μᵢ) Bₗ Πᵐₗ(μⱼ), where Bₗ contains the Greek
# coefficients and Π depends only on geometry. At fixed quadrature,
# dAᵐ/dp = Σₗ Πᵐₗ (dBₗ/dp) Πᵐₗ; see also Sanghavi, Davis & Eldering,
# JQSRT 133 (2014), (C.40). Thus this same accumulator handles B and dB.
# Code index l=1 denotes paper degree l=0. The final Fourier factor and
# backward (I,Q)↔(U,V) signs follow the package's Hovenier convention;
# they are identical for values and tangents (see conventions.md).
# Accumulate a whole Stokes block in registers instead of copying an Array
# slice for every l × angle pair × parameter. Round after each l as the
# reference accumulator does, including when μ has a wider floating type.
function _phase_column_tabulated!(Zplus, Zminus, mod, coefficients, m, Π_pair)
    Π_pair === nothing && return nothing
    Π, Πminus = Π_pair
    lmax = length(first(coefficients))
    Bs = [construct_B_matrix(mod, coefficients..., l) for l in (m+1):lmax]
    BT = eltype(Bs)
    nstokes = mod.n
    factor = m == 0 ? 1.0 : 2.0
    @inbounds for j in eachindex(Π[m+1]), i in eachindex(Π[m+1])
        plus, minus = zero(BT), zero(BT)
        for (k, l) in enumerate((m+1):lmax)
            plus = convert(BT, plus + Π[l][i] * Bs[k] * Π[l][j])
            minus = convert(BT, minus + Π[l][i] * Bs[k] * Πminus[l][j])
        end
        if nstokes == 1
            Zplus[i,j] = factor * plus
            Zminus[i,j] = factor * minus
        else
            for b in 1:nstokes, a in 1:nstokes
                row, col = (i-1)*nstokes+a, (j-1)*nstokes+b
                Zplus[row,col] = factor * plus[a,b]
                flip = (a <= 2 && b >= 3) || (a >= 3 && b <= 2)
                Zminus[row,col] = (flip ? -factor : factor) * minus[a,b]
            end
        end
    end
    return nothing
end

function _phase_spectral_column!(Zplus, Zminus, dZplus, dZminus,
                                  mod, values, tangents, m, Π_pair, s)
    coefficients = map(a -> view(a,:,s), values)
    _phase_column_tabulated!(view(Zplus,:,:,s), view(Zminus,:,:,s),
                             mod, coefficients, m, Π_pair)
    for p in axes(dZplus,1)
        coefficients = map(a -> view(a,p,:,s), tangents)
        _phase_column_tabulated!(view(dZplus,p,:,:,s), view(dZminus,p,:,:,s),
                                 mod, coefficients, m, Π_pair)
    end
    return nothing
end

function _compute_Z_moments_lin_tabulated(mod, μ, greek, lin_greek, m, arr_type)
    values = map(f -> getfield(greek,f), (:α,:β,:γ,:δ,:ϵ,:ζ))
    tangents = map(f -> getfield(lin_greek,f), (:α̇,:β̇,:γ̇,:δ̇,:ϵ̇,:ζ̇))
    β, dβ = values[2], tangents[2]
    spectral = ndims(β) == 2
    valid = (ndims(β) == 1 && ndims(dβ) == 2) ||
            (spectral && ndims(dβ) == 3 && size(dβ,3) == size(β,2))
    valid && size(dβ,2) == size(β,1) &&
        all(size(a) == size(β) for a in values) &&
        all(size(a) == size(dβ) for a in tangents) ||
        throw(DimensionMismatch("Greek values/tangents must be (l,)/(4,l) or (l,nSpec)/(4,l,nSpec)"))
    FT = eltype(β)
    ns = spectral ? size(β,2) : 1
    nb = length(μ) * mod.n
    Zplus, Zminus = zeros(FT,nb,nb,ns), zeros(FT,nb,nb,ns)
    dZplus, dZminus = zeros(FT,4,nb,nb,ns), zeros(FT,4,nb,nb,ns)
    # Π is independent of wavelength and microphysics. Only Bₗ and dBₗ
    # change across the spectral/parameter axes; build angular tables once.
    tables = ZMomentTables(μ, size(β,1))
    Π_pair = make_Π_lists(mod, tables, m)
    values = map(a -> reshape(collect(a),size(β,1),ns), values)
    tangents = map(a -> reshape(collect(a),4,size(β,1),ns), tangents)
    if Threads.nthreads() > 1 && ns >= 64
        Threads.@threads for s in 1:ns
            _phase_spectral_column!(Zplus,Zminus,dZplus,dZminus,
                                     mod,values,tangents,m,Π_pair,s)
        end
    else
        for s in 1:ns
            _phase_spectral_column!(Zplus,Zminus,dZplus,dZminus,
                                     mod,values,tangents,m,Π_pair,s)
        end
    end
    if spectral
        return arr_type(Zplus), arr_type(Zminus), arr_type(dZplus), arr_type(dZminus)
    end
    return arr_type(reshape(Zplus,nb,nb)), arr_type(reshape(Zminus,nb,nb)),
           arr_type(reshape(dZplus,4,nb,nb)), arr_type(reshape(dZminus,4,nb,nb))
end
