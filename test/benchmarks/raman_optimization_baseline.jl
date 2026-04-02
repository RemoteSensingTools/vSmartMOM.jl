#=
Raman Optimization Baseline — saves reference results for bit-exact comparison.
Output: test/benchmarks/raman_opttest_output/raman_reference.jld2
=#
using CUDA
device!(1)
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.SolarModel
using vSmartMOM.InelasticScattering
using Statistics, Interpolations, DelimitedFiles, JLD2

const ROOT = "/home/sanghavi/code/github/vSmartMOM.jl"
const OUTDIR = joinpath(ROOT, "test", "benchmarks", "raman_opttest_output")
mkpath(OUTDIR)

# Load 1-band config
parameters = parameters_from_yaml(joinpath(ROOT, "test/test_parameters/O2_parameters2_1band_opttest.yaml"))
FT = parameters.float_type
model = model_from_parameters(parameters)
model.obs_geom = CoreRT.ObsGeometry(parameters.sza, parameters.vza, parameters.vaz, parameters.obs_alt)
model.quad_points = CoreRT.rt_set_streams(parameters.quadrature_type, parameters.l_trunc,
    model.obs_geom, parameters.polarization_type, array_type(parameters.architecture))

# Raman setup
ν = model.params.spec_bands[1]
effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
n2, o2 = InelasticScattering.getRamanAtmoConstants(mean(ν), effT)

RS_type0 = InelasticScattering.noRS(fscattRayl=[FT(1)], ϖ_Cabannes=[FT(1)],
    bandSpecLim=[], iBand=[1], F₀=zeros(FT,1,1), SIF₀=zeros(FT,1,1))
RS_type1 = InelasticScattering.RRS(n2=n2, o2=o2,
    greek_raman=InelasticScattering.GreekCoefs([FT(1)],[FT(1)],[FT(1)],[FT(1)],[FT(1)],[FT(1)]),
    fscattRayl=[FT(1)], ϖ_Cabannes=[FT(1)], ϖ_λ₁λ₀=zeros(FT,1), i_λ₁λ₀=zeros(Int,1),
    Z⁻⁺_λ₁λ₀=zeros(FT,1,1), Z⁺⁺_λ₁λ₀=zeros(FT,1,1), i_ref=0, n_Raman=0,
    F₀=zeros(FT,1,1), SIF₀=zeros(FT,1,1))

CoreRT.getRamanSSProp!(RS_type1, 1e7/mean(ν), ν)
@info "n_Raman = $(RS_type1.n_Raman), nSpec = $(length(ν))"

# Solar spectrum
Tsolar = solar_transmission_from_file(joinpath(ROOT, "src/SolarModel/solar.out"))
Tsolar_interp = LinearInterpolation(Tsolar[4:end,1], Tsolar[4:end,2])
P = planck_spectrum_wn(5777.0, ν) * 2.1629e-05 * π
F₀ = Tsolar_interp.(ν) .* P
for rs in (RS_type0, RS_type1)
    rs.F₀ = zeros(model.params.polarization_type.n, length(P))
    rs.SIF₀ = zeros(model.params.polarization_type.n, length(ν))
    rs.F₀[1,:] = F₀
end

# Run RRS (single run — serves as both warmup and reference)
@info "Running RRS..."
t_rrs = @elapsed begin
    R1, T1, ieR1, ieT1 = CoreRT.rt_run_test(RS_type1, model, 1)
end
@info "RRS done: $(round(t_rrs, digits=1)) s"

# Run noRS
@info "Running noRS..."
t_nors = @elapsed begin
    R0, T0, _, _ = CoreRT.rt_run_test(RS_type0, model, 1)
end
@info "noRS done: $(round(t_nors, digits=1)) s"

# Save reference
ref_file = joinpath(OUTDIR, "raman_reference.jld2")
@save ref_file R1 T1 ieR1 ieT1 R0 T0 ν F₀ t_rrs t_nors
@info "Saved reference to $ref_file"
@info "RRS: $(round(t_rrs, digits=1))s, noRS: $(round(t_nors, digits=1))s, ratio: $(round(t_rrs/t_nors, digits=1))x"
