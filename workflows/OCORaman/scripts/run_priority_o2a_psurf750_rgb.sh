#!/usr/bin/env bash
set -euo pipefail

# Thin wrapper for the priority psurf=750 hPa RGB scenes.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PSURF="${PSURF:-750}"
exec "${SCRIPT_DIR}/run_priority_o2a_rgb.sh"
