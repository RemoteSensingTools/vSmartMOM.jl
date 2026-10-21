"""
    interpolate_phase_blocks(ν_spec, ν_nodes, values, AT)

Interpolate equal-shaped phase blocks directly into the requested backend.
Two nodes use linear interpolation; three use the same natural cubic spline
as `_interpolate_phase_nodes`. The spectral coordinate is always wavenumber.
Only the small nodal blocks and spectral grid are uploaded. There is no host
array per wavelength and no full spectral tensor copied from host to device.

The map is linear in the nodal values at fixed knots, so forward phase
matrices and their microphysical tangents use this identical operation.
This commutes with angular phase evaluation and preserves the upstream
truncated-Greek derivative, including the truncation-factor normalization.
"""
function interpolate_phase_blocks(ν_spec, ν_nodes, values, AT)
    length(ν_nodes) == length(values) || throw(DimensionMismatch("phase node/value count"))
    length(values) in (2,3) || throw(ArgumentError("phase interpolation requires 2 or 3 nodes"))
    p = sortperm(ν_nodes)
    knots = ν_nodes[p]
    all(diff(knots) .> 0) || throw(ArgumentError("phase interpolation knots must be distinct"))
    y = values[p]
    all(v->size(v)==size(first(y)),y) || throw(DimensionMismatch("phase node shapes"))
    FT = promote_type(eltype(ν_spec),eltype(knots),eltype(first(y)))
    grid = AT(FT.(ν_spec))
    nodes = map(v->AT(FT.(v)),y)
    out = similar(first(nodes),FT,size(first(y))...,length(ν_spec))
    flat = reshape(out,length(first(y)),length(ν_spec))
    backend = KernelAbstractions.get_backend(out)
    if length(nodes) == 2
        _interpolate_two_phase_nodes!(backend)(flat,grid,vec(nodes[1]),vec(nodes[2]),
            FT(knots[1]),FT(knots[2]);ndrange=size(flat))
    else
        h₀,h₁ = FT(knots[2]-knots[1]),FT(knots[3]-knots[2])
        # Natural boundary conditions: M₀=M₂=0. The sole interior second
        # derivative is M₁=3[(y₂-y₁)/h₁-(y₁-y₀)/h₀]/(h₀+h₁).
        # It is independent of evaluation wavelength and computed once.
        curvature = 3 .* ((nodes[3].-nodes[2])./h₁ .- (nodes[2].-nodes[1])./h₀) ./ (h₀+h₁)
        _interpolate_three_phase_nodes!(backend)(flat,grid,vec(nodes[1]),vec(nodes[2]),
            vec(nodes[3]),vec(curvature),FT(knots[1]),FT(knots[2]),FT(knots[3]);
            ndrange=size(flat))
    end
    return out
end

@kernel function _interpolate_two_phase_nodes!(out,@Const(grid),@Const(y₀),@Const(y₁),x₀,x₁)
    i,s = @index(Global,NTuple)
    @inbounds begin
        w = (grid[s]-x₀)/(x₁-x₀)
        out[i,s] = (one(w)-w)*y₀[i] + w*y₁[i]
    end
end

@kernel function _interpolate_three_phase_nodes!(out,@Const(grid),@Const(y₀),@Const(y₁),
        @Const(y₂),@Const(M₁),x₀,x₁,x₂)
    i,s = @index(Global,NTuple)
    @inbounds begin
        x = grid[s]
        if x <= x₁
            xa,xb,ya,yb,Ma,Mb = x₀,x₁,y₀[i],y₁[i],zero(M₁[i]),M₁[i]
        else
            xa,xb,ya,yb,Ma,Mb = x₁,x₂,y₁[i],y₂[i],M₁[i],zero(M₁[i])
        end
        h = xb-xa
        # Standard natural-cubic segment, matching _natural_cubic_three.
        out[i,s] = Ma*(xb-x)^3/(6h) + Mb*(x-xa)^3/(6h) +
            (ya-Ma*h^2/6)*(xb-x)/h + (yb-Mb*h^2/6)*(x-xa)/h
    end
end
