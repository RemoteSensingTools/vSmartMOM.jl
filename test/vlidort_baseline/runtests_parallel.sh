#!/usr/bin/env bash
# Parallel VLIDORT-baseline sweep: one Julia process per (quad, float, arch)
# spec. Each process caps BLAS to 8 threads; up to MAX_PARALLEL processes run
# concurrently. Walltime drops to (longest spec) instead of sum(specs).
#
# Usage:
#   cd test
#   ./vlidort_baseline/runtests_parallel.sh        # default: max 4 in parallel
#   MAX_PARALLEL=8 ./vlidort_baseline/runtests_parallel.sh
#
# Output: per-spec logs in /tmp/vlidort_baseline_v6_<spec>.log; combined
# summary printed at the end.

set -u
cd "$(dirname "$0")/.."   # → test/

MAX_PARALLEL="${MAX_PARALLEL:-4}"
LOGDIR="${LOGDIR:-/tmp/vlidort_baseline_parallel}"
mkdir -p "$LOGDIR"

SPECS=(
  "GaussQuadHemisphere/Float64/CPU"
  "GaussQuadHemisphere/Float64/GPU"
  "GaussQuadHemisphere/Float32/CPU"
  "GaussQuadHemisphere/Float32/GPU"
)

echo "Launching ${#SPECS[@]} specs, up to $MAX_PARALLEL in parallel"
echo "Logs: $LOGDIR/<spec>.log"
echo

start=$(date +%s)

# Use xargs to fan out with a parallelism cap. Each spec → one Julia invocation.
printf '%s\n' "${SPECS[@]}" | xargs -I{} -P "$MAX_PARALLEL" bash -c '
  spec="$1"
  safe_name=$(echo "$spec" | tr / _)
  log="'"$LOGDIR"'/${safe_name}.log"
  echo "[$(date +%H:%M:%S)] start  $spec"
  VLIDORT_BASELINE_ONLY_SPEC="$spec" VLIDORT_BASELINE_BLAS_THREADS=8 \
    julia --project=. -e "include(\"vlidort_baseline/runtests.jl\")" \
    > "$log" 2>&1
  rc=$?
  echo "[$(date +%H:%M:%S)] done   $spec (exit=$rc)"
' _ {}

end=$(date +%s)
echo
echo "Total walltime: $((end - start))s"
echo

# Combined per-spec summary table
echo "=== Per-spec summary ==="
for spec in "${SPECS[@]}"; do
  safe_name=$(echo "$spec" | tr / _)
  log="$LOGDIR/${safe_name}.log"
  if [[ ! -f "$log" ]]; then
    echo "  $spec: <no log>"
    continue
  fi
  grep -E "Case [ABC] —|Case [ABC] \[" "$log" | grep -E "max|TOA-up|BOA-dn" | head -8 \
    | sed "s/^/  /"
  echo "  ---"
done

# Final pass/fail aggregation
echo
echo "=== Aggregate Test Summary ==="
for spec in "${SPECS[@]}"; do
  safe_name=$(echo "$spec" | tr / _)
  log="$LOGDIR/${safe_name}.log"
  [[ ! -f "$log" ]] && continue
  line=$(grep -E "^VLIDORT baseline " "$log" | tail -1)
  echo "  $spec → ${line:-<no summary>}"
done
