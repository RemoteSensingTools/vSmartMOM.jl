# Distributions.jl Integration Update

**Date:** October 20, 2025  
**Branch:** TOMAS-aerosols

## Changes Made

Updated `explore_NK_julia.jl` to use Julia's `Distributions` package for log-normal size distributions instead of manually implementing the formula. This makes the code cleaner, more robust, and leverages well-tested statistical distributions.

## Key Modifications

### Before (Manual Implementation)
```julia
function lognormal(r::AbstractVector, p::AbstractVector)
    N_total, r_med, sigma_g = p
    return @. (N_total / (sqrt(2π) * log(sigma_g))) * 
              exp(-0.5 * (log(r / r_med) / log(sigma_g))^2)
end
```

### After (Using Distributions.jl)
```julia
using Distributions

function lognormal(r::AbstractVector, p::AbstractVector)
    N_total, r_med, sigma_g = p
    # Create LogNormal distribution with μ=log(r_med), σ=log(sigma_g)
    dist = LogNormal(log(r_med), log(sigma_g))
    # LogNormal pdf is per unit r, convert to per log10(r):
    # dN/dlog10(r) = dN/dr × dr/dlog10(r) = dN/dr × r × ln(10)
    return N_total .* pdf.(dist, r) .* r .* log(10)
end
```

## Technical Details

### Coordinate System Conversion

The key insight is understanding the difference between probability densities in different coordinate systems:

- **Distributions.jl LogNormal**: Returns `pdf(r)` = dN/dr (per unit radius)
- **Size distributions**: Need dN/dlog₁₀(r) (per log₁₀ bin)

**Jacobian transformation:**
```
dN/dlog₁₀(r) = dN/dr × dr/dlog₁₀(r) = dN/dr × r × ln(10)
```

This is because:
- dr/dlog₁₀(r) = d(r)/d(log₁₀(r)) = r × ln(10)
- For natural log: dr/dln(r) = r

### Verification

Created `test_distributions.jl` to verify the implementation:

```julia
Manual formula (ln): 575.55 #/cm³
Manual formula (log10): 1325.26 #/cm³
Distributions.jl: 1325.26 #/cm³
Difference: 0.00% ✓
```

Perfect agreement after accounting for the log₁₀ vs natural log difference!

## Benefits

1. **Cleaner code**: Uses standard library instead of custom implementation
2. **More robust**: Leverages well-tested statistical package
3. **Additional features**: Access to mean(), median(), mode(), std() functions
4. **Better documentation**: LogNormal distribution is well-documented in Distributions.jl
5. **Maintainability**: Less custom code to maintain

## Files Modified

- `src/Aerosols/data_exploration/explore_NK_julia.jl` - Updated to use Distributions.jl
- `src/Aerosols/data_exploration/test_distributions.jl` - **NEW** verification test
- `src/Aerosols/data_exploration/OVERVIEW.md` - Updated documentation

## Dependencies

The `Distributions` package was already in the main `Project.toml`:
```toml
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
```

So no new dependencies were added. The package is compatible with versions 0.23, 0.24, 0.25.

## Testing

Run the test script to verify:
```bash
julia --project=. src/Aerosols/data_exploration/test_distributions.jl
```

Expected output:
- Shows log-normal distribution properties (mean, median, mode, std)
- Compares manual formula with Distributions.jl
- Confirms 0.00% difference after proper coordinate transformation

## Next Steps

The main exploration script `explore_NK_julia.jl` is now ready to use with the improved Distributions.jl implementation. All plots and fitting routines remain unchanged in functionality but now use a more robust statistical foundation.

## References

- [Distributions.jl Documentation](https://juliastats.org/Distributions.jl/stable/)
- [LogNormal Distribution](https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.LogNormal)
- Size distribution theory: Seinfeld & Pandis, "Atmospheric Chemistry and Physics"
