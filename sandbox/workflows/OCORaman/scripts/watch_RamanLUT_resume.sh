#!/usr/bin/env bash
set -uo pipefail

REPO="${REPO:-/home/sanghavi/code/github/uni_vSmartMOM}"
PSURF="${PSURF:-750}"
CUDA_DEVICE="${CUDA_DEVICE:-1}"
OUT_NC="${OUT_NC:-/home/sanghavi/data/RamanSIFgrid/O2ABand/o2a_raman_lut_psurf${PSURF}.nc}"
INTERVAL_SECONDS="${INTERVAL_SECONDS:-600}"
LOGDIR="${LOGDIR:-$(dirname "$OUT_NC")}"
TAG="${TAG:-psurf${PSURF}}"
WATCH_LOG="${WATCH_LOG:-$LOGDIR/o2a_raman_lut_${TAG}_watchdog.log}"
PIDFILE="${PIDFILE:-$LOGDIR/o2a_raman_lut_${TAG}_watchdog.pid}"
LOCKFILE="${LOCKFILE:-$LOGDIR/o2a_raman_lut_${TAG}_watchdog.lock}"

mkdir -p "$LOGDIR"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')" "$*" | tee -a "$WATCH_LOG"
}

lut_process_running() {
    pgrep -u "$USER" -f 'julia .*createRamanLUT_O2A\.jl' >/dev/null 2>&1
}

lut_file_open() {
    fuser "$OUT_NC" >/dev/null 2>&1
}

lut_complete() {
    OUT_NC="$OUT_NC" julia --project="$REPO" -e '
        using NCDatasets
        const FILL_LIMIT = Float32(1e20)
        path = ENV["OUT_NC"]
        isfile(path) || exit(1)
        ds = NCDataset(path, "r")
        try
            alb = ds["albedo"][:]
            sza = ds["sza"][:]
            vars = (ds["stokes_rayleigh"], ds["stokes_cabannes"], ds["stokes_rrs"])
            for ia in eachindex(alb), isza in eachindex(sza)
                for A in vars
                    vals = Array(A[1, 1, ia, isza, :, :])
                    physical = isfinite.(vals) .& (abs.(vals) .< FILL_LIMIT)
                    (all(physical) && maximum(abs.(vals)) > 0f0) || exit(1)
                end
            end
            exit(0)
        finally
            close(ds)
        end
    ' >/dev/null 2>&1
}

start_resume() {
    local run_log="$LOGDIR/o2a_raman_lut_${TAG}_resume_$(date '+%Y%m%d_%H%M%S').log"
    log "Starting RESUME=1 PSURFS=$PSURF CUDA_DEVICE=$CUDA_DEVICE OUT_NC=$OUT_NC"
    setsid bash -lc "cd '$REPO' && exec env RESUME=1 CUDA_DEVICE='$CUDA_DEVICE' PSURFS='$PSURF' OUT_NC='$OUT_NC' julia --project=. workflows/OCORaman/scripts/createRamanLUT_O2A.jl" \
        > "$run_log" 2>&1 < /dev/null &
    local pid="$!"
    log "Started PID=$pid log=$run_log"
}

main_loop() {
    echo "$$" > "$PIDFILE"
    log "Watchdog started for $TAG; interval=${INTERVAL_SECONDS}s"
    log "PIDFILE=$PIDFILE LOCKFILE=$LOCKFILE"

    while true; do
        if lut_process_running; then
            log "LUT process is running; no action"
        elif lut_file_open; then
            log "NetCDF is open by another process; no action"
        elif lut_complete; then
            log "LUT appears complete; watchdog exiting"
            rm -f "$PIDFILE"
            exit 0
        else
            log "No active LUT process and LUT incomplete; restarting"
            start_resume
        fi
        sleep "$INTERVAL_SECONDS"
    done
}

(
    flock -n 9 || {
        printf '[%s] Another watchdog already holds %s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')" "$LOCKFILE" >&2
        exit 1
    }
    main_loop
) 9>"$LOCKFILE"
