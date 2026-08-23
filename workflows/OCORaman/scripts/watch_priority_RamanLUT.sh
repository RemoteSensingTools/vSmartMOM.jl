#!/usr/bin/env bash
set -uo pipefail

# Conservative watchdog + timing logger for priority O2A Raman LUT chunks.
#
# It monitors the psurf-specific priority wrapper:
#   run_priority_o2a_psurf1000_rgb.sh  or  run_priority_o2a_psurf750_rgb.sh
#
# Behavior:
#   - records periodic progress/GPU/log snapshots
#   - rewrites a parsed RT-scene timing CSV from chunk logs
#   - restarts the priority wrapper only when no matching driver/JULIA process
#     is running and the priority run is not complete
#   - does not kill or restart a live process just because a long RT step is quiet
#
# Typical:
#   PSURF=750 CUDA_HOME=/usr/local/cuda CUDA_DEVICE=1 \
#   nohup workflows/OCORaman/scripts/watch_priority_RamanLUT.sh \
#     > /home/sanghavi/data/RamanSIFgrid/O2ABand/o2a_raman_lut_priority_psurf0750_rgb/logs/watchdog.out 2>&1 &
#
#   PSURF=1000 CUDA_HOME=/usr/local/cuda CUDA_DEVICE=1 \
#   nohup workflows/OCORaman/scripts/watch_priority_RamanLUT.sh \
#     > /home/sanghavi/data/RamanSIFgrid/O2ABand/o2a_raman_lut_priority_psurf1000_rgb/logs/watchdog.out 2>&1 &

REPO="${REPO:-/home/sanghavi/code/github/uni_vSmartMOM}"
PSURF="${PSURF:-750}"
CUDA_DEVICE="${CUDA_DEVICE:-1}"
CUDA_HOME="${CUDA_HOME:-/usr/local/cuda}"
INTERVAL_SECONDS="${INTERVAL_SECONDS:-300}"
STALE_WARN_SECONDS="${STALE_WARN_SECONDS:-1800}"
JULIA_FLAGS="${JULIA_FLAGS:---pkgimages=no}"
JULIA_BIN="${JULIA_BIN:-julia}"
PRIORITY_ALBEDO_IDXS="${PRIORITY_ALBEDO_IDXS:-1,7,21}"

PSURF_TAG="$(printf '%04.0f' "${PSURF}")"
PSURF_WRAPPER_TAG="$(printf '%.0f' "${PSURF}")"
OUTDIR="${OUTDIR:-/home/sanghavi/data/RamanSIFgrid/O2ABand/o2a_raman_lut_priority_psurf${PSURF_TAG}_rgb}"
LOGDIR="${LOGDIR:-${OUTDIR}/logs}"
WATCH_LOG="${WATCH_LOG:-${LOGDIR}/priority_psurf${PSURF_TAG}_watchdog.log}"
SNAPSHOT_CSV="${SNAPSHOT_CSV:-${LOGDIR}/priority_psurf${PSURF_TAG}_timing_snapshots.csv}"
EVENT_CSV="${EVENT_CSV:-${LOGDIR}/priority_psurf${PSURF_TAG}_rt_scene_timings.csv}"
PIDFILE="${PIDFILE:-${LOGDIR}/priority_psurf${PSURF_TAG}_watchdog.pid}"
LOCKFILE="${LOCKFILE:-${LOGDIR}/priority_psurf${PSURF_TAG}_watchdog.lock}"
NOHUP_DRIVER_LOG="${NOHUP_DRIVER_LOG:-${LOGDIR}/nohup_driver.log}"
WRAPPER="${WRAPPER:-workflows/OCORaman/scripts/run_priority_o2a_psurf${PSURF_WRAPPER_TAG}_rgb.sh}"

mkdir -p "${LOGDIR}" "${OUTDIR}"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')" "$*" | tee -a "${WATCH_LOG}"
}

init_csvs() {
    if [[ ! -f "${SNAPSHOT_CSV}" ]]; then
        printf 'epoch,iso,host,psurf,pid,pid_elapsed,driver_running,gpu_util_percent,gpu_mem_mib,gpu_power_w,log_age_s,log_size_bytes,nc_count,completed_chunk_count,rt_done_count,current_chunk,current_sza,current_albedo,current_band,note\n' > "${SNAPSHOT_CSV}"
    fi
}

driver_running() {
    local pid
    for pid in $(pgrep -u "${USER}" -f 'bash .*run_priority_o2a_rgb\.sh|bash .*run_priority_o2a_psurf[0-9]+_rgb\.sh' 2>/dev/null); do
        [[ -r "/proc/${pid}/environ" ]] || continue
        tr '\0' '\n' < "/proc/${pid}/environ" | grep -qx "PSURF=${PSURF}" && return 0
    done
    return 1
}

julia_pid() {
    local pid
    for pid in $(pgrep -u "${USER}" -f 'julia .*chunked_createRamanLUT_O2A\.jl' 2>/dev/null); do
        [[ -r "/proc/${pid}/environ" ]] || continue
        if tr '\0' '\n' < "/proc/${pid}/environ" | grep -qx "OUTDIR=${OUTDIR}"; then
            printf '%s\n' "${pid}"
            return 0
        fi
    done
    return 1
}

run_complete() {
    local sza_idx missing
    for sza_idx in $(seq 1 14); do
        missing="$(
            env \
                PSURF="${PSURF}" \
                SZA_IDX="${sza_idx}" \
                ALBEDO_IDXS="${PRIORITY_ALBEDO_IDXS}" \
                OUTDIR="${OUTDIR}" \
                "${JULIA_BIN}" ${JULIA_FLAGS} --project=. workflows/OCORaman/scripts/priority_missing_RamanLUT_O2A.jl \
                2>/dev/null | tail -n 1
        )"
        [[ -n "${missing}" ]] && return 1
    done
    return 0
}

latest_log_file() {
    find "${LOGDIR}" -maxdepth 1 -type f -name "priority_psurf${PSURF_TAG}_sza*_alb*.log" \
        -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -n 1 | cut -d' ' -f2-
}

current_from_logs() {
    local log_file="$1"
    local chunk="" sza="" albedo="" band=""
    if [[ -n "${log_file}" && -f "${log_file}" ]]; then
        chunk="$(rg -o 'priority_psurf[0-9]+_sza[0-9]+_alb[0-9_,-]+' "${log_file}" 2>/dev/null | tail -n 1)"
        sza="$(rg -n 'sza = ' "${log_file}" 2>/dev/null | tail -n 1 | sed -E 's/.*sza = ([^[:space:]]+).*/\1/')"
        albedo="$(rg -n 'albedo = ' "${log_file}" 2>/dev/null | tail -n 1 | sed -E 's/.*albedo = ([^[:space:]]+).*/\1/')"
        band="$(rg -n 'ib = ' "${log_file}" 2>/dev/null | tail -n 1 | sed -E 's/.*ib = ([^[:space:]]+).*/\1/')"
    fi
    printf '%s,%s,%s,%s\n' "${chunk}" "${sza}" "${albedo}" "${band}"
}

gpu_snapshot() {
    nvidia-smi --id="${CUDA_DEVICE}" \
        --query-gpu=utilization.gpu,memory.used,power.draw \
        --format=csv,noheader,nounits 2>/dev/null |
        head -n 1 |
        awk -F',' '{gsub(/^ +| +$/, "", $1); gsub(/^ +| +$/, "", $2); gsub(/^ +| +$/, "", $3); printf "%s,%s,%s", $1, $2, $3}'
}

write_timing_events() {
    {
        printf 'host,psurf,log_file,sza,albedo,sif_on,band,seconds\n'
        awk -v host="$(hostname)" -v psurf_watch="${PSURF}" '
            /RT scene done/ {
                active=1; psurf=""; sza=""; albedo=""; sif=""; band=""; seconds="";
                next;
            }
            active && /psurf =/ { psurf=$NF; next; }
            active && /sza =/ { sza=$NF; next; }
            active && /albedo =/ { albedo=$NF; next; }
            active && /sif_on =/ { sif=$NF; next; }
            active && /ib =/ { band=$NF; next; }
            active && /seconds =/ {
                seconds=$NF;
                gsub(/[^0-9.]/, "", seconds);
                if (seconds != "") {
                    printf "%s,%s,%s,%s,%s,%s,%s,%s\n", host, psurf_watch, FILENAME, sza, albedo, sif, band, seconds;
                }
                active=0;
            }
        ' "${LOGDIR}"/priority_psurf"${PSURF_TAG}"_sza*_alb*.log 2>/dev/null
    } > "${EVENT_CSV}.tmp"
    mv "${EVENT_CSV}.tmp" "${EVENT_CSV}"
}

append_snapshot() {
    local pid="$1" note="$2"
    local epoch iso host elapsed driver_flag gpu log_file log_mtime log_age log_size nc_count chunk_count done_count current
    epoch="$(date +%s)"
    iso="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    host="$(hostname)"
    elapsed=""
    if [[ -n "${pid}" ]]; then
        elapsed="$(ps -p "${pid}" -o etime= 2>/dev/null | awk '{$1=$1; print}')"
    fi
    if driver_running; then driver_flag=1; else driver_flag=0; fi
    gpu="$(gpu_snapshot)"
    [[ -n "${gpu}" ]] || gpu=",,"
    log_file="$(latest_log_file)"
    log_mtime=0
    log_size=0
    if [[ -n "${log_file}" && -f "${log_file}" ]]; then
        log_mtime="$(stat -c '%Y' "${log_file}")"
        log_size="$(stat -c '%s' "${log_file}")"
    fi
    log_age=$(( epoch - log_mtime ))
    nc_count="$(find "${OUTDIR}" -maxdepth 1 -type f -name "priority_psurf${PSURF_TAG}_sza*_alb*.nc" 2>/dev/null | wc -l | awk '{print $1}')"
    chunk_count="$(grep -c "Finished priority O2A Raman LUT run for psurf=${PSURF}" "${LOGDIR}/driver.log" 2>/dev/null || true)"
    done_count="$(rg -c 'RT scene done' "${LOGDIR}"/priority_psurf"${PSURF_TAG}"_sza*_alb*.log 2>/dev/null | awk -F: '{s += $2} END {print s+0}')"
    current="$(current_from_logs "${log_file}")"
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "${epoch}" "${iso}" "${host}" "${PSURF}" "${pid}" "${elapsed}" "${driver_flag}" \
        "${gpu}" "${log_age}" "${log_size}" "${nc_count}" "${chunk_count}" "${done_count}" \
        "${current}" "${note}" >> "${SNAPSHOT_CSV}"
}

start_priority_run() {
    log "Starting priority wrapper PSURF=${PSURF} CUDA_DEVICE=${CUDA_DEVICE} JULIA_FLAGS=${JULIA_FLAGS}"
    setsid bash -lc "cd '${REPO}' && exec env CUDA_HOME='${CUDA_HOME}' CUDA_DEVICE='${CUDA_DEVICE}' JULIA_FLAGS='${JULIA_FLAGS}' '${WRAPPER}'" \
        >> "${NOHUP_DRIVER_LOG}" 2>&1 < /dev/null &
    log "Started priority wrapper PID=$!; output=${NOHUP_DRIVER_LOG}"
}

main_loop() {
    echo "$$" > "${PIDFILE}"
    init_csvs
    log "Priority watchdog started for psurf=${PSURF}; interval=${INTERVAL_SECONDS}s stale_warn=${STALE_WARN_SECONDS}s"
    log "OUTDIR=${OUTDIR}"
    log "SNAPSHOT_CSV=${SNAPSHOT_CSV}"
    log "EVENT_CSV=${EVENT_CSV}"

    while true; do
        local pid="" note="" latest="" age=0
        pid="$(julia_pid || true)"
        latest="$(latest_log_file)"
        if [[ -n "${latest}" && -f "${latest}" ]]; then
            age=$(( $(date +%s) - $(stat -c '%Y' "${latest}") ))
        fi

        write_timing_events

        if [[ -n "${pid}" ]]; then
            note="running"
            if (( age > STALE_WARN_SECONDS )); then
                note="running_log_stale_${age}s"
                log "Warning: process PID=${pid} is running, but latest log has not changed for ${age}s"
            fi
            append_snapshot "${pid}" "${note}"
            log "Running PID=${pid}; note=${note}; latest_log_age=${age}s"
        elif driver_running; then
            note="driver_running_no_julia"
            append_snapshot "" "${note}"
            log "Driver wrapper is running but no matching Julia process is visible; waiting"
        elif run_complete; then
            note="complete"
            append_snapshot "" "${note}"
            log "Priority run appears complete; watchdog exiting"
            rm -f "${PIDFILE}"
            exit 0
        else
            note="restart"
            append_snapshot "" "${note}"
            log "No priority process found and run is incomplete; restarting"
            start_priority_run
        fi

        sleep "${INTERVAL_SECONDS}"
    done
}

(
    flock -n 9 || {
        printf '[%s] Another priority watchdog already holds %s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${LOCKFILE}" >&2
        exit 1
    }
    main_loop
) 9>"${LOCKFILE}"
