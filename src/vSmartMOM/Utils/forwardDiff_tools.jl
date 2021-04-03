

# Define batched matrix multiply for GPU and Duals:
function batched_mul(A::CuArray{ForwardDiff.Dual{T,V,N},3}, B::CuArray{ForwardDiff.Dual{T,V,N},3}) where {T,V,N}
    # Extract values:
    Av = ForwardDiff.value.(A)
    Bv = ForwardDiff.value.(B)
    # Use strided batch for A*B (defined as gemm_strided_batched):
    Cv = Av ⊠ Bv
    # Compute derivatives ∂(AB)/∂x = A * ∂B/∂x + ∂A/∂x * B;
    dABdx = [Av ⊠ ForwardDiff.partials.(B,i) + ForwardDiff.partials.(A,i) ⊠ Bv for i=1:N];
    dABdx = ForwardDiff.Partials.(tuple.(dABdx...));
    return eltype(A).(Cv,dABdx);
end


# Overload of batch_inv! for Dual numbers
function batch_inv!(X::CuArray{ForwardDiff.Dual{T,V,N},3}, A::CuArray{ForwardDiff.Dual{T,V,N},3}) where {T,V,N}
    Atemp = ForwardDiff.value.(A)
    invA = similar(Atemp);
    # Set invA=A⁻¹
    batch_inv!(invA,Atemp)
    # Compute derivatives ∂A⁻¹/∂x = -A⁻¹ * ∂A/∂x * A⁻¹; using NNlib batched matrix multiply
    # @show typeof(ForwardDiff.partials.(A,1))
    -invA ⊠ ForwardDiff.partials.(A,1) ⊠ invA
    dAdx = [-invA ⊠ ForwardDiff.partials.(A,i) ⊠ invA for i=1:N];
    #println("Pack into tuples again")
    dAdx = ForwardDiff.Partials.(tuple.(dAdx...));
    X .= eltype(X).(invA,dAdx);

end