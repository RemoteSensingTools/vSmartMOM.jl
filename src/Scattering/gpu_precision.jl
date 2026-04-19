#=
GPU-compatible precision primitives for Mie scattering computations.

Three precision tiers:
  Tier 1 -- DoubleSingle / ComplexDS: emulated FP64 via pairs of Float32
            (~48 mantissa bits). Used for the Dn downward recursion where
            Float32 alone loses significant digits for |y| > 50.
  Tier 2 -- NeumaierAccum: compensated summation in native Float32.
            Eliminates systematic drift in long running sums (S1/S2
            amplitudes, cross-section dot products, Greek-coefficient
            quadrature).
  Tier 3 -- Native Float32 (no helpers needed).

All functions are @inline and allocation-free for GPU kernel compatibility.
=#

# ============================================================================
# Precision policy dispatch -- selects Tier 1 strategy per GPU capability
# ============================================================================

"""Abstract precision policy for Mie GPU kernels."""
abstract type MiePrecisionPolicy end

"""Use native Float64 for the Dn recursion (A100, V100 -- full FP64 throughput)."""
struct NativeFloat64 <: MiePrecisionPolicy end

"""Use DoubleSingle emulation for the Dn recursion (L40S, consumer GPUs -- FP32 only)."""
struct DSEmulated    <: MiePrecisionPolicy end

# ============================================================================
# DoubleSingle{T} -- error-free arithmetic building blocks
# ============================================================================

"""
    DoubleSingle{T}

Unevaluated sum of two floating-point numbers `hi + lo` where `|lo| <= ulp(hi)/2`.
Provides ~2p mantissa bits (48 bits for T=Float32).

All operations are branchless and use only `+`, `-`, `*`, `fma`.
"""
struct DoubleSingle{T <: AbstractFloat}
    hi::T
    lo::T
end

@inline DoubleSingle(x::T) where {T} = DoubleSingle{T}(x, zero(T))
@inline DoubleSingle{T}(x::T) where {T} = DoubleSingle{T}(x, zero(T))

# Convert from higher precision for testing
@inline function DoubleSingle{T}(x::Float64) where {T <: AbstractFloat}
    hi = T(x)
    lo = T(x - Float64(hi))
    DoubleSingle{T}(hi, lo)
end

@inline Base.convert(::Type{T}, ds::DoubleSingle{T}) where {T} = ds.hi + ds.lo
@inline Base.convert(::Type{Float64}, ds::DoubleSingle{Float32}) = Float64(ds.hi) + Float64(ds.lo)

# --- Error-free transformations ---

"""
    TwoSum(a, b) -> DoubleSingle{T}

Error-free addition: returns (s, e) such that a + b = s + e exactly.
Knuth's algorithm (6 FLOPs).
"""
@inline function TwoSum(a::T, b::T) where {T <: AbstractFloat}
    s = a + b
    v = s - a
    e = (a - (s - v)) + (b - v)
    DoubleSingle{T}(s, e)
end

"""
    TwoProd(a, b) -> DoubleSingle{T}

Error-free multiplication using hardware `fma`: returns (p, e) such that
a * b = p + e exactly (2 FLOPs with fma).
"""
@inline function TwoProd(a::T, b::T) where {T <: AbstractFloat}
    p = a * b
    e = fma(a, b, -p)
    DoubleSingle{T}(p, e)
end

# --- DoubleSingle arithmetic ---

"""Add two DoubleSingle values."""
@inline function ds_add(a::DoubleSingle{T}, b::DoubleSingle{T}) where {T}
    s = TwoSum(a.hi, b.hi)
    # Accumulate low-order terms into the error
    e = s.lo + a.lo + b.lo
    # Re-normalize
    TwoSum(s.hi, e)
end

"""Subtract two DoubleSingle values."""
@inline function ds_sub(a::DoubleSingle{T}, b::DoubleSingle{T}) where {T}
    ds_add(a, DoubleSingle{T}(-b.hi, -b.lo))
end

"""Multiply two DoubleSingle values."""
@inline function ds_mul(a::DoubleSingle{T}, b::DoubleSingle{T}) where {T}
    p = TwoProd(a.hi, b.hi)
    # Cross terms (only need hi parts for ~2p accuracy)
    e = p.lo + a.hi * b.lo + a.lo * b.hi
    TwoSum(p.hi, e)
end

"""Divide two DoubleSingle values: a / b."""
@inline function ds_div(a::DoubleSingle{T}, b::DoubleSingle{T}) where {T}
    # First approximation
    q1 = a.hi / b.hi
    # Compute remainder: r = a - q1 * b
    r = ds_sub(a, ds_mul(DoubleSingle{T}(q1), b))
    # Correction
    q2 = r.hi / b.hi
    TwoSum(q1, q2)
end

"""Reciprocal of a DoubleSingle value: 1 / a."""
@inline function ds_inv(a::DoubleSingle{T}) where {T}
    ds_div(DoubleSingle{T}(one(T), zero(T)), a)
end

# --- Convenience: DS ↔ scalar ---

@inline ds_add(a::DoubleSingle{T}, b::T) where {T} = ds_add(a, DoubleSingle{T}(b, zero(T)))
@inline ds_add(a::T, b::DoubleSingle{T}) where {T} = ds_add(DoubleSingle{T}(a, zero(T)), b)
@inline ds_sub(a::DoubleSingle{T}, b::T) where {T} = ds_sub(a, DoubleSingle{T}(b, zero(T)))
@inline ds_sub(a::T, b::DoubleSingle{T}) where {T} = ds_sub(DoubleSingle{T}(a, zero(T)), b)
@inline ds_mul(a::DoubleSingle{T}, b::T) where {T} = ds_mul(a, DoubleSingle{T}(b, zero(T)))
@inline ds_mul(a::T, b::DoubleSingle{T}) where {T} = ds_mul(DoubleSingle{T}(a, zero(T)), b)
@inline ds_div(a::DoubleSingle{T}, b::T) where {T} = ds_div(a, DoubleSingle{T}(b, zero(T)))
@inline ds_div(a::T, b::DoubleSingle{T}) where {T} = ds_div(DoubleSingle{T}(a, zero(T)), b)

# ============================================================================
# ComplexDS{T} -- complex arithmetic in DoubleSingle precision
# ============================================================================

"""
    ComplexDS{T}

Complex number with DoubleSingle real and imaginary parts (4 floats of type T).
Used for the Dn recursion where complex division has catastrophic cancellation.
"""
struct ComplexDS{T <: AbstractFloat}
    re::DoubleSingle{T}
    im::DoubleSingle{T}
end

@inline ComplexDS(re::DoubleSingle{T}, im::DoubleSingle{T}) where {T} = ComplexDS{T}(re, im)

"""Construct from a single real scalar (imaginary part zero)."""
@inline cds_real(x::T) where {T <: AbstractFloat} = ComplexDS{T}(DoubleSingle{T}(x), DoubleSingle{T}(zero(T)))

"""Construct from real and imaginary scalars."""
@inline cds_complex(re::T, im::T) where {T <: AbstractFloat} = ComplexDS{T}(DoubleSingle{T}(re), DoubleSingle{T}(im))

# Convert from Complex{Float64} for testing
@inline function ComplexDS{T}(z::Complex{Float64}) where {T <: AbstractFloat}
    ComplexDS{T}(DoubleSingle{T}(T(real(z)), T(real(z) - Float64(T(real(z))))),
                 DoubleSingle{T}(T(imag(z)), T(imag(z) - Float64(T(imag(z))))))
end

# Convert back to Complex{Float64}
@inline function Base.convert(::Type{Complex{Float64}}, z::ComplexDS{Float32})
    Complex{Float64}(convert(Float64, z.re), convert(Float64, z.im))
end

# Convert to native complex
@inline function to_complex(z::ComplexDS{T}) where {T}
    Complex{T}(z.re.hi + z.re.lo, z.im.hi + z.im.lo)
end

"""Add two ComplexDS values."""
@inline function cds_add(a::ComplexDS{T}, b::ComplexDS{T}) where {T}
    ComplexDS{T}(ds_add(a.re, b.re), ds_add(a.im, b.im))
end

"""Subtract two ComplexDS values."""
@inline function cds_sub(a::ComplexDS{T}, b::ComplexDS{T}) where {T}
    ComplexDS{T}(ds_sub(a.re, b.re), ds_sub(a.im, b.im))
end

"""Multiply two ComplexDS values: (a.re + i*a.im)(b.re + i*b.im)."""
@inline function cds_mul(a::ComplexDS{T}, b::ComplexDS{T}) where {T}
    # (ar*br - ai*bi) + i*(ar*bi + ai*br)
    re = ds_sub(ds_mul(a.re, b.re), ds_mul(a.im, b.im))
    im = ds_add(ds_mul(a.re, b.im), ds_mul(a.im, b.re))
    ComplexDS{T}(re, im)
end

"""Divide two ComplexDS values: a / b = a * conj(b) / |b|^2."""
@inline function cds_div(a::ComplexDS{T}, b::ComplexDS{T}) where {T}
    # |b|^2 = b.re^2 + b.im^2
    denom = ds_add(ds_mul(b.re, b.re), ds_mul(b.im, b.im))
    # a * conj(b) = (ar*br + ai*bi) + i*(ai*br - ar*bi)
    num_re = ds_add(ds_mul(a.re, b.re), ds_mul(a.im, b.im))
    num_im = ds_sub(ds_mul(a.im, b.re), ds_mul(a.re, b.im))
    ComplexDS{T}(ds_div(num_re, denom), ds_div(num_im, denom))
end

"""Reciprocal of ComplexDS: 1 / z."""
@inline function cds_inv(z::ComplexDS{T}) where {T}
    one_ds = ComplexDS{T}(DoubleSingle{T}(one(T)), DoubleSingle{T}(zero(T)))
    cds_div(one_ds, z)
end

# --- Convenience: ComplexDS ↔ scalar ---

@inline function cds_add(a::ComplexDS{T}, b::T) where {T}
    ComplexDS{T}(ds_add(a.re, b), a.im)
end
@inline function cds_add(a::T, b::ComplexDS{T}) where {T}
    ComplexDS{T}(ds_add(a, b.re), b.im)
end

"""Multiply ComplexDS by real DoubleSingle scalar."""
@inline function cds_mul_real(a::ComplexDS{T}, s::DoubleSingle{T}) where {T}
    ComplexDS{T}(ds_mul(a.re, s), ds_mul(a.im, s))
end

"""Construct ComplexDS from a real DoubleSingle (imaginary = 0)."""
@inline function cds_from_real(s::DoubleSingle{T}) where {T}
    ComplexDS{T}(s, DoubleSingle{T}(zero(T)))
end

# ============================================================================
# NeumaierAccum{T} -- compensated summation (Tier 2)
# ============================================================================

"""
    NeumaierAccum{T}

Neumaier's improved Kahan summation. Maintains a running sum `s` and a
compensation variable `c`. Unlike plain Kahan, this handles the case where
`|addend| > |s|` correctly.

Cost: ~7 FLOPs per addition (vs 1 for naive sum).
"""
struct NeumaierAccum{T <: AbstractFloat}
    s::T   # running sum
    c::T   # compensation for lost low-order bits
end

@inline NeumaierAccum{T}() where {T} = NeumaierAccum{T}(zero(T), zero(T))
@inline NeumaierAccum(s::T) where {T} = NeumaierAccum{T}(s, zero(T))

"""Add a value to the Neumaier accumulator."""
@inline function neumaier_add(acc::NeumaierAccum{T}, x::T) where {T}
    t = acc.s + x
    # If |s| >= |x|, low-order bits of x are lost; otherwise low-order bits of s
    c = if abs(acc.s) >= abs(x)
        acc.c + ((acc.s - t) + x)
    else
        acc.c + ((x - t) + acc.s)
    end
    NeumaierAccum{T}(t, c)
end

"""Extract the compensated sum."""
@inline neumaier_sum(acc::NeumaierAccum{T}) where {T} = acc.s + acc.c

# --- Complex Neumaier for S1/S2 accumulations ---

"""
    ComplexNeumaier{T}

Pair of Neumaier accumulators for complex summation.
"""
struct ComplexNeumaier{T <: AbstractFloat}
    re::NeumaierAccum{T}
    im::NeumaierAccum{T}
end

@inline ComplexNeumaier{T}() where {T} = ComplexNeumaier{T}(NeumaierAccum{T}(), NeumaierAccum{T}())

"""Add a complex value to the complex Neumaier accumulator."""
@inline function cneumaier_add(acc::ComplexNeumaier{T}, re::T, im::T) where {T}
    ComplexNeumaier{T}(neumaier_add(acc.re, re), neumaier_add(acc.im, im))
end

"""Extract the compensated complex sum."""
@inline function cneumaier_sum(acc::ComplexNeumaier{T}) where {T}
    Complex{T}(neumaier_sum(acc.re), neumaier_sum(acc.im))
end
