"""Run the Julia benchmark and sample its process's reserved GPU memory.

Usage: python3 monitor_gpu.py JULIA_EXECUTABLE GPU_INDEX LABEL [POLARIZATION]
Run from test/. Results go to AUDIT_OUTPUT (an existing directory).
This samples nvidia-smi process memory, including the CUDA memory pool; it
does not measure exact live array bytes or guarantee capture of short peaks.
"""
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

julia, gpu, label = sys.argv[1:4]
pol = sys.argv[4] if len(sys.argv) > 4 else "Stokes_I()"
out = Path(os.environ["AUDIT_OUTPUT"])
log = out / (label + ".log")
env = dict(os.environ, AUDIT_BACKEND="cuda", AUDIT_NSPEC="10000",
           AUDIT_POL=pol, CUDA_VISIBLE_DEVICES=gpu,
           JULIA_NUM_THREADS="4", OPENBLAS_NUM_THREADS="1")
script = Path(__file__).resolve().with_name("benchmark.jl")
with log.open("w") as stream:
    process = subprocess.Popen([julia, "--startup-file=no", "--project=.", str(script)],
                               stdout=stream, stderr=subprocess.STDOUT, env=env)
    peaks = {}
    try:
        while process.poll() is None:
            phases = re.findall(r"^BENCH_PHASE (.+)$", log.read_text(), re.M)
            phase = phases[-1] if phases else "construction/compilation"
            result = subprocess.run(
                ["nvidia-smi", "--id=" + gpu,
                 "--query-compute-apps=pid,used_gpu_memory", "--format=csv,noheader,nounits"],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                universal_newlines=True, check=True)
            for line in result.stdout.splitlines():
                pid, memory = (part.strip() for part in line.split(",", 1))
                if pid == str(process.pid):
                    peaks[phase] = max(peaks.get(phase, 0), int(memory))
            time.sleep(0.25)
    finally:
        if process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
    summary = {"exit_code": process.returncode, "gpu_index": gpu,
               "metric": "sampled process GPU memory in MiB, including reserved pool",
               "sample_interval_seconds": 0.25, "peak_mib_by_phase": peaks}
    (out / (label + "-memory.json")).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))
    sys.exit(process.returncode)
