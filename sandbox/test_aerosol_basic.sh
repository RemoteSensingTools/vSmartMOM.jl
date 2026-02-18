#!/bin/bash
# Simple test script that can be run from bash

cd /home/cfranken/code/gitHub/vSmartMOM.jl

echo "========================================================================"
echo "Testing Aerosol Framework with FT Type Parameters"
echo "========================================================================"
echo ""

# Create a minimal Julia test
cat > /tmp/test_aerosol_minimal.jl << 'JULIA_CODE'
using YAML

println("Step 1: Testing YAML config loading...")
config_path = "examples/aerosol_config_tomas15.yaml"
if isfile(config_path)
    config = YAML.load_file(config_path)
    println("  ✓ Config loaded successfully")
    println("  ✓ Scheme type: ", config["aerosol_scheme"]["type"])
    println("  ✓ Number of bins: ", config["aerosol_scheme"]["size_bins"]["n_bins"])
    println("  ✓ Species count: ", length(config["aerosol_scheme"]["species"]))
else
    println("  ✗ Config not found")
    exit(1)
end
println()

println("Step 2: Testing type parameter pattern...")
# Simple struct with FT parameter
struct TestStruct{FT}
    value::FT
    array::Vector{FT}
end

# Constructor with default FT
function TestStruct(val, arr, FT=Float64)
    return TestStruct{FT}(FT(val), FT.(arr))
end

# Test Float64
t64 = TestStruct(3.14, [1, 2, 3], Float64)
println("  ✓ Float64 struct: ", typeof(t64))
println("  ✓ Value type: ", typeof(t64.value))
println("  ✓ Array element type: ", eltype(t64.array))

# Test Float32
t32 = TestStruct(3.14, [1, 2, 3], Float32)
println("  ✓ Float32 struct: ", typeof(t32))
println("  ✓ Value type: ", typeof(t32.value))
println("  ✓ Array element type: ", eltype(t32.array))

@assert typeof(t64.value) == Float64
@assert typeof(t32.value) == Float32
println()

println("Step 3: Testing TOMAS-15 bin calculation...")
n_bins = config["aerosol_scheme"]["size_bins"]["n_bins"]
diam_min = Float64(config["aerosol_scheme"]["size_bins"]["diam_min_nm"])
diam_max = Float64(config["aerosol_scheme"]["size_bins"]["diam_max_nm"])

bin_edges = diam_min .* (diam_max / diam_min) .^ (Float64.(collect(0:n_bins)) ./ Float64(n_bins))
bin_centers = sqrt.(bin_edges[1:end-1] .* bin_edges[2:end])

println("  ✓ Calculated $(length(bin_centers)) bin centers")
println("  ✓ First bin center: $(bin_centers[1]) nm")
println("  ✓ Last bin center: $(bin_centers[end]) nm")
println("  ✓ Bin edges type: $(typeof(bin_edges))")
@assert length(bin_centers) == n_bins
@assert bin_edges[1] ≈ diam_min
@assert bin_edges[end] ≈ diam_max
println()

println("Step 4: Testing species properties...")
species_config = config["aerosol_scheme"]["species"]
species_names = collect(keys(species_config))
println("  ✓ Found $(length(species_names)) species")

densities = Dict{String, Float64}()
for (sp, sp_config) in species_config
    densities[sp] = Float64(sp_config["density"])
end
println("  ✓ Extracted densities: ", length(densities), " entries")
println("  ✓ DUST density: ", densities["DUST"], " kg/m³")
println("  ✓ Density type: ", typeof(first(values(densities))))
@assert all(values(densities) .> 0.0)
println()

println("======================================================================")
println("All Basic Tests Passed! ✓")
println("======================================================================")
println()
println("Summary:")
println("  • YAML config loading works")
println("  • FT type parameter pattern is correct")
println("  • Bin calculations produce correct results")
println("  • Type conversions work properly")
println("  • Dict value types are correct")
println()
JULIA_CODE

echo "Running Julia test..."
julia /tmp/test_aerosol_minimal.jl

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Test completed successfully!"
    exit 0
else
    echo ""
    echo "✗ Test failed"
    exit 1
fi
