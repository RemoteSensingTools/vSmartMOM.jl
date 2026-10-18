using Test, vSmartMOM, YAML, Logging
using vSmartMOM.CoreRT

@testset "Shared parsed HITRAN data" begin
    clear_spectroscopy_cache!()
    tasks = [Threads.@spawn CoreRT._hitran_lines("O2",Float64) for _ in 1:4]
    cached = fetch.(tasks)
    @test all(x -> x === first(cached), cached)
    @test length(CoreRT._HITRAN_LINES) == 1
    fresh = CoreRT.AtmosphericAbsorption.load_lines(
        CoreRT.AtmosphericAbsorption.HitranPort(artifact("O2")); FT=Float64)
    @test first(cached).ν0 == fresh.ν0
    @test first(cached).S == fresh.S
    @test CoreRT._hitran_lines("O2",Float32) !== first(cached)
    @test length(CoreRT._HITRAN_LINES) == 2
    with_logger(NullLogger()) do
        cfg = YAML.load_file("test_parameters/JacobianTestFast.yaml")
        delete!(cfg,"scattering")
        params = read_parameters(cfg)
        params.spec_bands[1] = collect(range(first(params.spec_bands[1]),last(params.spec_bands[1]);length=5))
        forward = model_from_parameters(params)
        linearized,_ = model_from_parameters(LinMode(),params)
        @test forward.τ_abs[1] ≈ linearized.τ_abs[1] rtol=1e-12
        @test CoreRT._hitran_lines("O2",Float64) === first(cached)
        @test first(cached).S == fresh.S
        context = BatchContext(params)
        @test only(only(context.absorption_models)).lines === first(cached)
        @test context.model.τ_abs[1] ≈ forward.τ_abs[1] rtol=1e-12
    end
    clear_spectroscopy_cache!()
    @test isempty(CoreRT._HITRAN_LINES)
    @test first(cached).S == fresh.S # clearing ownership does not mutate users
end

# --- Downloaded-edition cache identity and coverage (offline) ---------------
#
# Every test below runs in a temporary scratch directory with the hitran.org
# download replaced by a local fixture server; nothing touches the network.
# The fixtures are synthetic 160-column HITRAN records: they exercise cache
# selection, not spectroscopic accuracy.

using Printf, Scratch

function hitran_record(mol, iso, ν; S=1.0e-22)
    record = lpad(mol, 2) * string(iso) * @sprintf("%12.6f", ν) * @sprintf("%10.3E", S) *
             @sprintf("%10.3E", 1.0) * ".0700" * "0.090" * @sprintf("%10.4f", 100.0) *
             "0.75" * @sprintf("%8.6f", 0.0) * " "^60 * "000000" * " "^12 * " " *
             "    1.0" * "    1.0"
    @assert ncodeunits(record) == 160
    return record
end
hitran_file(records...) = join(records, "\n") * "\n"

const CO2_IDS = Int[7, 8, 9, 10, 11, 12, 13, 14, 121, 15, 120, 122]
const CO2_NARROW = hitran_file(hitran_record(2, 1, 6200.0), hitran_record(2, 2, 6201.0))
const CO2_FULL = hitran_file(hitran_record(2, 1, 2349.0), hitran_record(2, 1, 6200.0),
                             hitran_record(2, 2, 6201.0))

"""Write a cache entry and its metadata as `fetch_hitran_by_ids` does."""
function seed_cache!(name, edition, content; numin=0, numax=150000,
                     iso_ids=CO2_IDS, drop=(), overrides...)
    par = vSmartMOM._hitran_par_path(name, edition)
    mkpath(dirname(par))
    write(par, content)
    fields = Dict{Symbol,Any}(:molecule => name, :edition => edition,
        :download_date => "2026-10-17T00:00:00", :numin => numin, :numax => numax,
        :iso_ids => join(iso_ids, ","), :source_url => "fixture",
        :sha256 => bytes2hex(vSmartMOM.sha256(content)),
        :file_size_bytes => ncodeunits(content))
    merge!(fields, Dict(overrides))
    order = (:molecule, :edition, :download_date, :numin, :numax, :iso_ids,
             :source_url, :sha256, :file_size_bytes)
    vSmartMOM._write_meta_toml(vSmartMOM._hitran_meta_path(name, edition);
        (k => fields[k] for k in order if !(k in drop))...)
    return par
end

"""Offline stand-in for hitran.org: records each URL, then runs `respond`."""
struct FixtureServer
    urls::Vector{String}
    respond::Function
    lock::ReentrantLock
end
FixtureServer(urls, respond) = FixtureServer(urls, respond, ReentrantLock())
function (server::FixtureServer)(url, path)
    lock(() -> push!(server.urls, url), server.lock)
    server.respond(path)
    return path
end
serve(content) = FixtureServer(String[], path -> write(path, content))
failing(err) = FixtureServer(String[], path -> throw(err))
const OFFLINE = ErrorException("network access is disabled in tests")

function with_offline_hitran(f)
    previous_download = vSmartMOM._HITRAN_DOWNLOAD[]
    previous_edition = get_hitran_edition()
    vSmartMOM._HITRAN_DOWNLOAD[] = failing(OFFLINE)
    try
        mktempdir() do root
            Scratch.with_scratch_directory(root) do
                with_logger(f, NullLogger())
            end
        end
    finally
        vSmartMOM._HITRAN_DOWNLOAD[] = previous_download
        set_hitran_edition!(previous_edition)
        clear_spectroscopy_cache!()
    end
end

function with_server(f, server)
    previous = vSmartMOM._HITRAN_DOWNLOAD[]
    vSmartMOM._HITRAN_DOWNLOAD[] = server
    try
        return f()
    finally
        vSmartMOM._HITRAN_DOWNLOAD[] = previous
    end
end

cache_files(edition) = sort(readdir(dirname(vSmartMOM._hitran_par_path("CO2", edition))))

@testset "HITRAN download cache identity and coverage" begin with_offline_hitran() do

@testset "requested isotopologues must be included (reproduction)" begin
    # Metadata for global ID 7 alone must not satisfy a request for [7, 8].
    seed_cache!("CO2", "iso", hitran_file(hitran_record(2, 1, 6200.0));
                numin=6000, numax=6400, iso_ids=[7])
    server = serve(CO2_NARROW)
    path = with_server(server) do
        fetch_hitran_by_ids("CO2", [7, 8]; numin=6000, numax=6400, edition="iso")
    end
    @test length(server.urls) == 1
    @test occursin("iso_ids_list=7,8&", only(server.urls))
    @test read(path, String) == CO2_NARROW
    @test vSmartMOM._read_meta_toml(vSmartMOM._hitran_meta_path("CO2", "iso"))["iso_ids"] == "7,8"
    # A file holding more isotopologues than requested is reused.
    server = serve(CO2_NARROW)
    with_server(server) do
        fetch_hitran_by_ids("CO2", [8]; numin=6100, numax=6300, edition="iso")
    end
    @test isempty(server.urls)
end

@testset "spectral coverage" begin
    seed_cache!("CO2", "band", CO2_NARROW; numin=6000, numax=6400)
    server = serve(CO2_NARROW)
    with_server(() -> fetch_hitran("CO2"; numin=6100, numax=6300, edition="band"), server)
    @test isempty(server.urls)
    server = serve(CO2_NARROW)
    with_server(() -> fetch_hitran("CO2"; numin=5900, numax=6300, edition="band"), server)
    @test length(server.urls) == 1
end

@testset "artifact() rejects a narrower scratch file (reproduction)" begin
    set_hitran_edition!("narrow")
    seed_cache!("CO2", "narrow", CO2_NARROW; numin=6000, numax=6400)
    @test !hitran_is_cached("CO2", "narrow")
    server = serve(CO2_FULL)
    path = with_server(() -> artifact("CO2"), server)
    @test length(server.urls) == 1
    @test occursin("&numin=0&numax=150000", only(server.urls))
    @test read(path, String) == CO2_FULL
    @test hitran_is_cached("CO2", "narrow")
    # A verified full cache is used without downloading (unchanged default).
    server = serve(CO2_NARROW)
    @test with_server(() -> artifact("CO2"), server) == path
    @test isempty(server.urls)
    @test read(path, String) == CO2_FULL
end

@testset "artifact() fails clearly instead of reusing unverified data" begin
    set_hitran_edition!("offline")
    par = seed_cache!("CO2", "offline", CO2_NARROW; numin=6000, numax=6400)
    err = try
        artifact("CO2"); nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("does not cover [0, 150000]", err.msg)
    @test occursin("set_hitran_edition!(\"artifact\")", err.msg)
    @test read(par, String) == CO2_NARROW
end

@testset "missing or inconsistent provenance is not verified (reproduction)" begin
    status(edition) = vSmartMOM._hitran_cache_status("CO2", edition, CO2_IDS, 0, 150000)
    cases = [
        ("nometa",    "missing",           () -> rm(vSmartMOM._hitran_meta_path("CO2", "nometa"))),
        ("sha",       "SHA-256",           () -> seed_cache!("CO2", "sha", CO2_FULL; sha256="0"^64)),
        ("size",      "file size",         () -> seed_cache!("CO2", "size", CO2_FULL; file_size_bytes=1)),
        ("nosha",     "sha256",            () -> seed_cache!("CO2", "nosha", CO2_FULL; drop=(:sha256,))),
        ("noiso",     "iso_ids",           () -> seed_cache!("CO2", "noiso", CO2_FULL; drop=(:iso_ids,))),
        ("badiso",    "iso_ids",           () -> seed_cache!("CO2", "badiso", CO2_FULL; iso_ids=["7", "x"])),
        ("molecule",  "molecule",          () -> seed_cache!("CO2", "molecule", CO2_FULL; molecule="H2O")),
        ("edition",   "edition",           () -> seed_cache!("CO2", "edition", CO2_FULL; edition="HITRAN2016")),
        ("nobounds",  "numin",             () -> seed_cache!("CO2", "nobounds", CO2_FULL; drop=(:numin,))),
        ("bigiso",    "iso_ids",           () -> (meta = vSmartMOM._hitran_meta_path("CO2", "bigiso");
                                                  write(meta, replace(read(meta, String), r"iso_ids = .*" => "iso_ids = 1e100")))),
        # A crash between the two renames leaves new data beside old metadata.
        ("crash",     "SHA-256",           () -> write(vSmartMOM._hitran_par_path("CO2", "crash"),
                                                       replace(CO2_FULL, "6201.000000" => "6202.000000"))),
    ]
    for (edition, _, corrupt!) in cases
        # fetch_hitran re-downloads rather than reusing.
        seed_cache!("CO2", edition, CO2_FULL)
        corrupt!()
        @test !hitran_is_cached("CO2", edition)
        server = serve(CO2_FULL)
        with_server(() -> fetch_hitran("CO2"; edition=edition), server)
        @test length(server.urls) == 1
        @test hitran_is_cached("CO2", edition)
        # artifact() does the same.
        seed_cache!("CO2", edition, CO2_FULL)
        corrupt!()
        set_hitran_edition!(edition)
        server = serve(CO2_FULL)
        with_server(() -> artifact("CO2"), server)
        @test length(server.urls) == 1
    end
    # The shared predicate names the failed check.
    for (edition, reason, corrupt!) in cases
        seed_cache!("CO2", edition, CO2_FULL)
        corrupt!()
        verified, why = status(edition)
        @test !verified
        @test occursin(reason, why)
    end
end

@testset "content rewrites are detected after a verified use" begin
    par = seed_cache!("CO2", "rewrite", CO2_FULL)
    @test hitran_is_cached("CO2", "rewrite")
    write(par, replace(CO2_FULL, "6201.000000" => "6202.000000"))  # same size
    @test !hitran_is_cached("CO2", "rewrite")
    # Also after the hash of an older file has been remembered.
    par = seed_cache!("CO2", "aged", CO2_FULL)
    sleep(2.5)
    @test hitran_is_cached("CO2", "aged")
    @test haskey(vSmartMOM._HITRAN_HASHES, par)
    write(par, replace(CO2_FULL, "6201.000000" => "6202.000000"))
    @test !hitran_is_cached("CO2", "aged")
end

@testset "an unreadable cache is reported, not thrown" begin
    par = seed_cache!("CO2", "unreadable", CO2_FULL)
    chmod(par, 0o000)
    try
        # Permissions cannot block reading on Windows or for root.
        blocked = try
            read(par); false
        catch
            true
        end
        if blocked
            verified, why = vSmartMOM._hitran_cache_status("CO2", "unreadable", CO2_IDS, 0, 150000)
            @test !verified
            @test occursin("could not be checked", why)
            @test !hitran_is_cached("CO2", "unreadable")
        else
            @test_skip "file permissions do not block reading here"
        end
    finally
        chmod(par, 0o644)
    end
end

@testset "failed downloads preserve the previous valid cache (reproduction)" begin
    edition = "atomic"
    par = seed_cache!("CO2", edition, CO2_FULL)
    meta = vSmartMOM._hitran_meta_path("CO2", edition)
    old_meta = read(meta, String)
    responses = [
        ("transport error after partial data", path -> (write(path, CO2_FULL[1:100]); error("connection reset"))),
        ("empty body",                         path -> write(path, "")),
        ("truncated record",                   path -> write(path, CO2_FULL[1:end-20])),
        ("malformed numeric field",            path -> write(path, hitran_file(replace(hitran_record(2, 1, 6200.0), "6200.000000" => "6200.0000xx")))),
        ("unrequested isotopologue",           path -> write(path, hitran_file(hitran_record(1, 1, 6200.0)))),
        ("line outside requested band",        path -> write(path, hitran_file(hitran_record(2, 1, 7000.0)))),
        ("record longer than 160 columns",     path -> write(path, hitran_file(hitran_record(2, 1, 6200.0) * " "))),
        ("records joined without a newline",   path -> write(path, hitran_file(hitran_record(2, 1, 6200.0) * hitran_record(2, 2, 6201.0)))),
    ]
    for (label, respond) in responses
        server = FixtureServer(String[], respond)
        @test_throws Exception with_server(server) do
            fetch_hitran("CO2"; numin=6000, numax=6400, edition=edition, force=true)
        end
        @test length(server.urls) == 1
        @test read(par, String) == CO2_FULL
        @test read(meta, String) == old_meta
        @test cache_files(edition) == ["CO2.meta.toml", "CO2.par"]
        @test hitran_is_cached("CO2", edition)
    end
    # A successful replacement leaves exactly the new verified pair.
    server = serve(CO2_NARROW)
    with_server(server) do
        fetch_hitran("CO2"; numin=6000, numax=6400, edition=edition, force=true)
    end
    @test read(par, String) == CO2_NARROW
    @test cache_files(edition) == ["CO2.meta.toml", "CO2.par"]
    @test vSmartMOM._hitran_cache_status("CO2", edition, CO2_IDS, 6000, 6400)[1]
end

@testset "transport errors fail even with a 2xx status" begin
    check = vSmartMOM._check_hitran_response
    url = "https://hitran.org/lbl/api"
    ok = vSmartMOM.Downloads.Response("https", url, 200, "HTTP/1.1 200 OK", Pair{String,String}[])
    @test check(ok) === nothing
    partial = vSmartMOM.Downloads.RequestError(url, 18, "transfer closed with outstanding read data remaining", ok)
    @test_throws ErrorException check(partial)
    @test_throws ErrorException check(vSmartMOM.Downloads.Response("https", url, 500, "", Pair{String,String}[]))
    limited = vSmartMOM.Downloads.Response("https", url, 403, "", Pair{String,String}[])
    @test_throws r"rate limit" check(limited)
end

@testset "failed replacement restores the previous pair" begin
    edition = "replace"
    par = seed_cache!("CO2", edition, CO2_FULL)
    meta = vSmartMOM._hitran_meta_path("CO2", edition)
    old_meta = read(meta, String)
    for failing_target in (par, meta)
        tmp_par, tmp_meta = tempname(dirname(par)), tempname(dirname(par))
        write(tmp_par, CO2_NARROW)
        write(tmp_meta, "replacement metadata\n")
        rename = (src, dst) -> dst == failing_target ? error("injected rename failure") :
                                                       vSmartMOM._atomic_rename(src, dst)
        @test_throws ErrorException vSmartMOM._replace_hitran_cache!(
            tmp_par, par, tmp_meta, meta; rename=rename)
        @test read(par, String) == CO2_FULL
        @test read(meta, String) == old_meta
        @test hitran_is_cached("CO2", edition)
        rm(tmp_par; force=true); rm(tmp_meta; force=true)
        @test cache_files(edition) == ["CO2.meta.toml", "CO2.par"]
    end
    # Once the new pair is in place, failing to delete the backup is only a warning.
    tmp_par, tmp_meta = tempname(dirname(par)), tempname(dirname(par))
    write(tmp_par, CO2_NARROW)
    write(tmp_meta, "replacement metadata\n")
    remove = (path; force) -> error("injected cleanup failure")
    @test vSmartMOM._replace_hitran_cache!(tmp_par, par, tmp_meta, meta; remove=remove) == par
    @test read(par, String) == CO2_NARROW
    @test read(meta, String) == "replacement metadata\n"
end

@testset "concurrent fetches publish a consistent pair" begin
    # Both downloads finish before either is published.
    arrived = Threads.Atomic{Int}(0)
    server = FixtureServer(String[], path -> begin
        first_arrival = Threads.atomic_add!(arrived, 1) == 0
        write(path, first_arrival ? CO2_NARROW : hitran_file(hitran_record(2, 1, 6300.0)))
        deadline = time() + 60
        while arrived[] < 2 && time() < deadline
            yield()
        end
    end)
    with_server(server) do
        tasks = [Threads.@spawn fetch_hitran("CO2"; numin=6000, numax=6400,
                                             edition="concurrent", force=true) for _ in 1:2]
        foreach(wait, tasks)
    end
    @test length(server.urls) == 2
    @test vSmartMOM._hitran_cache_status("CO2", "concurrent", CO2_IDS, 6000, 6400)[1]
    @test cache_files("concurrent") == ["CO2.meta.toml", "CO2.par"]
end

@testset "managed replacement invalidates parsed spectroscopy (reproduction)" begin
    set_hitran_edition!("parsed")
    seed_cache!("CO2", "parsed", CO2_FULL)
    clear_spectroscopy_cache!()
    @test CoreRT._hitran_lines("CO2", Float64).ν0 == [2349.0, 6200.0, 6201.0]
    replacement = hitran_file(hitran_record(2, 1, 4000.0))
    with_server(() -> fetch_hitran("CO2"; edition="parsed", force=true), serve(replacement))
    @test CoreRT._hitran_lines("CO2", Float64).ν0 == [4000.0]
end

end end
