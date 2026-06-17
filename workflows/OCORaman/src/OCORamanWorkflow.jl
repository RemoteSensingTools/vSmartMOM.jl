module OCORamanWorkflow

using Dates
using JLD2
using JSON
using Printf
using Statistics

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering

export workflow_root,
       repo_root,
       output_root,
       figures_root,
       ensure_dir,
       git_sha,
       git_status_short,
       run_id,
       write_json,
       write_manifest,
       latest_jld2,
       first_lastdim_series,
       summarize_grid_file,
       print_grid_summary,
       raman_air_summary,
       print_raman_air_summary

workflow_root() = normpath(joinpath(@__DIR__, ".."))
repo_root() = normpath(joinpath(workflow_root(), "..", ".."))
output_root() = get(ENV, "OCORAMAN_OUTPUT_ROOT", joinpath(workflow_root(), "output"))
figures_root() = get(ENV, "OCORAMAN_FIGURES_ROOT", joinpath(workflow_root(), "figures"))

ensure_dir(path::AbstractString) = (mkpath(path); path)

function git_sha(root::AbstractString=repo_root())
    try
        return chomp(read(`git -C $root rev-parse HEAD`, String))
    catch
        return "unknown"
    end
end

function git_status_short(root::AbstractString=repo_root())
    try
        return chomp(read(`git -C $root status --short`, String))
    catch
        return ""
    end
end

function run_id(prefix::AbstractString="ocoraman")
    stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    sha = git_sha()
    short = sha == "unknown" ? "unknown" : first(sha, min(8, lastindex(sha)))
    return "$(prefix)_$(stamp)_$(short)"
end

function write_json(path::AbstractString, data)
    ensure_dir(dirname(path))
    open(path, "w") do io
        JSON.print(io, data, 2)
        println(io)
    end
    return path
end

function write_manifest(outdir::AbstractString; script="", yaml="", extra=Dict{String,Any}())
    ensure_dir(outdir)
    manifest = Dict{String,Any}(
        "created_at" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
        "workflow_root" => workflow_root(),
        "repo_root" => repo_root(),
        "git_sha" => git_sha(),
        "git_status_short" => git_status_short(),
        "script" => script,
        "yaml" => yaml,
        "extra" => extra,
    )
    return write_json(joinpath(outdir, "run_manifest.json"), manifest)
end

function latest_jld2(dir::AbstractString)
    files = filter(path -> endswith(lowercase(path), ".jld2"), readdir(dir; join=true))
    isempty(files) && error("No .jld2 files found in $(dir)")
    return sort(files; by=mtime)[end]
end

function first_lastdim_series(A)
    ndims(A) == 0 && return [A]
    idx = ntuple(d -> d == ndims(A) ? Colon() : 1, ndims(A))
    return vec(A[idx...])
end

function _maybe_read(f, key::AbstractString)
    key in keys(f) || return nothing
    return f[key]
end

function summarize_grid_file(path::AbstractString)
    jldopen(path, "r") do f
        nu = _maybe_read(f, "ν")
        r_rrs = _maybe_read(f, "R_rrs")
        r_nors = _maybe_read(f, "R_nors")
        s_rrs = r_rrs === nothing ? nothing : first_lastdim_series(r_rrs)
        s_nors = r_nors === nothing ? nothing : first_lastdim_series(r_nors)
        delta = (s_rrs === nothing || s_nors === nothing) ? nothing : s_rrs .- s_nors
        return (
            path = path,
            keys = collect(keys(f)),
            n_spectral = nu === nothing ? missing : length(nu),
            wavelength_nm_minmax = nu === nothing ? missing : extrema(1.0e7 ./ nu),
            R_rrs_size = r_rrs === nothing ? missing : size(r_rrs),
            R_nors_size = r_nors === nothing ? missing : size(r_nors),
            R_delta_minmax = delta === nothing ? missing : extrema(delta),
            R_delta_mean = delta === nothing ? missing : mean(delta),
        )
    end
end

function print_grid_summary(path::AbstractString)
    s = summarize_grid_file(path)
    println("Grid file: ", s.path)
    println("  n_spectral: ", s.n_spectral)
    println("  wavelength_nm_minmax: ", s.wavelength_nm_minmax)
    println("  R_rrs_size: ", s.R_rrs_size)
    println("  R_nors_size: ", s.R_nors_size)
    println("  R_rrs - R_nors min/max: ", s.R_delta_minmax)
    println("  R_rrs - R_nors mean: ", s.R_delta_mean)
    return s
end

function raman_air_summary(; wavelength_nm::Real=760.0,
                             temperature_K::Real=300.0,
                             vmr_n2::Real=0.8,
                             vmr_o2::Real=0.2)
    nu_eff = 1.0e7 / Float64(wavelength_nm)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(
        nu_eff, Float64(temperature_K), Float64(vmr_n2), Float64(vmr_o2))

    gamma_rayleigh, sigma_rayleigh =
        InelasticScattering.compute_γ_air_Rayleigh!(Float64(wavelength_nm), n2, o2)
    gamma_cabannes, omega_cabannes =
        InelasticScattering.compute_γ_air_Cabannes!(Float64(wavelength_nm), n2, o2)

    omega_n2, gamma_cab_n2, gamma_rayl_n2 =
        InelasticScattering.compute_γ_mol_Cabannes!(Float64(wavelength_nm), n2)
    omega_o2, gamma_cab_o2, gamma_rayl_o2 =
        InelasticScattering.compute_γ_mol_Cabannes!(Float64(wavelength_nm), o2)

    return (
        wavelength_nm = Float64(wavelength_nm),
        temperature_K = Float64(temperature_K),
        vmr_n2 = Float64(vmr_n2),
        vmr_o2 = Float64(vmr_o2),
        gamma_air_rayleigh = gamma_rayleigh,
        gamma_air_cabannes = gamma_cabannes,
        omega_air_cabannes = omega_cabannes,
        sigma_air_rayleigh_cm2_per_molecule = sigma_rayleigh,
        gamma_n2_rayleigh = gamma_rayl_n2,
        gamma_n2_cabannes = gamma_cab_n2,
        omega_n2_cabannes = omega_n2,
        gamma_o2_rayleigh = gamma_rayl_o2,
        gamma_o2_cabannes = gamma_cab_o2,
        omega_o2_cabannes = omega_o2,
    )
end

function print_raman_air_summary(; kwargs...)
    s = raman_air_summary(; kwargs...)
    @printf("λ = %.4f nm, T = %.2f K\n", s.wavelength_nm, s.temperature_K)
    @printf("  air γ Rayleigh: %.10f\n", s.gamma_air_rayleigh)
    @printf("  air γ Cabannes: %.10f\n", s.gamma_air_cabannes)
    @printf("  air ϖ Cabannes: %.10f\n", s.omega_air_cabannes)
    @printf("  air σ Rayleigh: %.8e cm^2 molecule^-1\n",
            s.sigma_air_rayleigh_cm2_per_molecule)
    @printf("  N2 γ Rayleigh/Cabannes/ϖ: %.10f / %.10f / %.10f\n",
            s.gamma_n2_rayleigh, s.gamma_n2_cabannes, s.omega_n2_cabannes)
    @printf("  O2 γ Rayleigh/Cabannes/ϖ: %.10f / %.10f / %.10f\n",
            s.gamma_o2_rayleigh, s.gamma_o2_cabannes, s.omega_o2_cabannes)
    return s
end

end
