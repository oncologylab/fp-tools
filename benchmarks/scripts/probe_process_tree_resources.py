#!/usr/bin/env python3
"""Qualify Linux delegated-cgroup accounting using synthetic workloads only.

No scientific data or model is read. This is an instrumentation test, not a
timing result for fp-tools. A later frozen benchmark may reuse measure_tree.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import uuid


def kernel_fields(path):
    return dict(line.split() for line in Path(path).read_text().splitlines())


def current_cgroup():
    lines = Path("/proc/self/cgroup").read_text().splitlines()
    unified = [line[3:] for line in lines if line.startswith("0::/")]
    if len(unified) != 1:
        raise RuntimeError("one unified cgroup v2 membership is required")
    return Path("/sys/fs/cgroup") / unified[0].lstrip("/")


def counter(path):
    value = Path(path).read_text().strip()
    if not value.isdecimal():
        raise ValueError(f"invalid kernel counter: {Path(path).name}")
    return int(value)


def measure_tree(command, cwd, environment, root, output):
    """Measure a fresh workload subtree of our dedicated delegated service.

    The supervisor moves outside the measured subtree before spawning. The
    child enters the subtree before exec; its descendants inherit membership.
    Empty cgroups persist until counters are read, including memory.peak.
    """
    root, output = Path(root), Path(output)
    if root != current_cgroup() or not root.name.startswith("fp-tools-resource-probe-"):
        raise RuntimeError("refusing to modify another cgroup or a non-probe service")
    if not command or not all(isinstance(arg, str) for arg in command):
        raise ValueError("command must be a nonempty argument array")
    supervisor, workload = root / "supervisor", root / "workload"
    supervisor.mkdir()
    (supervisor / "cgroup.procs").write_text(str(os.getpid()))
    (root / "cgroup.subtree_control").write_text("+cpu +memory")
    workload.mkdir()
    if counter(workload / "memory.peak") != 0:
        raise RuntimeError("workload cgroup is not fresh")
    if int(kernel_fields(workload / "cpu.stat")["usage_usec"]) != 0:
        raise RuntimeError("workload CPU counter is not fresh")

    def enter_workload():
        (workload / "cgroup.procs").write_text(str(os.getpid()))

    # This supervisor is single-threaded. preexec_fn only moves its forked
    # child into the owned cgroup before exec, including all command imports.
    started = time.monotonic()
    with (output / "workload.log").open("xb") as log:
        process = subprocess.Popen(command, cwd=cwd, env=environment,
                                   stdout=log, stderr=subprocess.STDOUT,
                                   preexec_fn=enter_workload)
        status = process.wait()
    elapsed = time.monotonic() - started
    populated = kernel_fields(workload / "cgroup.events")["populated"] != "0"
    cpu = kernel_fields(workload / "cpu.stat")
    peak = counter(workload / "memory.peak")
    result = {
        "schema": "fp-tools-cgroup-measurement-v1", "command": command,
        "returncode": status, "descendants_remaining": populated,
        "runtime_seconds": elapsed, "cpu_seconds": int(cpu["usage_usec"]) / 1e6,
        "user_cpu_seconds": int(cpu["user_usec"]) / 1e6,
        "system_cpu_seconds": int(cpu["system_usec"]) / 1e6,
        "peak_memory_bytes": peak, "peak_memory_mb": peak / 1e6,
        "cpu_scope": "process_tree", "memory_scope": "absolute_process_tree_peak",
        "memory_metric": "Linux_cgroup_v2_memory.peak",
        "memory_units": "decimal_MB",
        "memory_interpretation": "Kernel-accounted memory including charged file cache and kernel memory; not RSS. The same counter and cache policy must be used for both methods.",
        "measurement_complete": status == 0 and not populated,
        "supervisor_excluded": True,
    }
    (output / "measurement.json").write_text(json.dumps(result, indent=2) + "\n")
    if populated:
        # Only our newly created workload subtree is affected.
        (workload / "cgroup.kill").write_text("1")
        raise RuntimeError("workload left descendants alive; measurement rejected")
    workload.rmdir()
    if status:
        raise RuntimeError(f"workload failed with exit {status}; measurement retained")
    return result


def probe_command(mode):
    child = ("import time; x=bytearray(32*1024*1024); start=time.process_time();\n"
             "while time.process_time()-start<0.2: pass\n"
             "time.sleep(0.6)")
    if mode == "serial":
        body = "[subprocess.run([sys.executable,'-c',code],check=True) for _ in range(2)]"
    elif mode == "parallel":
        body = ("children=[subprocess.Popen([sys.executable,'-c',code]) for _ in range(2)]; "
                "statuses=[p.wait() for p in children]; assert statuses==[0,0]")
    elif mode == "failure":
        return [sys.executable, "-c", "raise SystemExit(7)"]
    elif mode == "orphan":
        return [sys.executable, "-c",
                "import subprocess,sys; subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)'])"]
    else:
        raise ValueError("unknown synthetic probe")
    return [sys.executable, "-c", f"import subprocess,sys; code={child!r}; {body}"]


def worker(output, mode):
    output = Path(output)
    root = current_cgroup()
    if not root.name.startswith("fp-tools-resource-probe-") or not root.name.endswith(".service"):
        raise RuntimeError("worker must run in its own fp-tools probe service")
    # Constrain the entire tree, not just numerical-library threads.
    cpus = sorted(os.sched_getaffinity(0))[:2]
    if len(cpus) != 2:
        raise RuntimeError("two process-visible CPUs required for concurrency probe")
    os.sched_setaffinity(0, cpus)
    environment = dict(os.environ)
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        environment[name] = "1"
    return measure_tree(probe_command(mode), output, environment, root, output)


def run_probe(outdir):
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=False)
    records = {}
    for mode in ("serial", "parallel", "failure", "orphan"):
        target = outdir / mode
        target.mkdir()
        unit = "fp-tools-resource-probe-" + uuid.uuid4().hex
        command = ["systemd-run", "--user", f"--unit={unit}", "--wait",
                   "--property=Type=exec", "--property=Delegate=yes",
                   "--property=MemoryAccounting=yes", "--property=CPUAccounting=yes",
                   "--property=RuntimeMaxSec=90", "--property=KillMode=control-group",
                   sys.executable, str(Path(__file__).resolve()),
                   "--worker", mode, "--outdir", str(target)]
        with (target / "launcher.log").open("x") as log:
            try:
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=110)
            finally:
                subprocess.run(["systemctl", "--user", "stop", unit],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
                subprocess.run(["systemctl", "--user", "reset-failed", unit],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        record = json.loads((target / "measurement.json").read_text())
        records[mode] = record
        if mode == "failure":
            if result.returncode == 0 or record["returncode"] != 7 or record["measurement_complete"]:
                raise RuntimeError("failed workload was not rejected")
        elif mode == "orphan":
            if (result.returncode == 0 or record["returncode"] != 0
                    or not record["descendants_remaining"] or record["measurement_complete"]):
                raise RuntimeError("orphan workload was not rejected")
        elif result.returncode or not record["measurement_complete"]:
            raise RuntimeError(f"{mode} probe failed; inspect retained logs")
    validate_probes(records)
    summary = {"schema": "fp-tools-process-tree-instrumentation-probe-v1",
               "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "kernel": os.uname().release,
               "systemd": subprocess.check_output(["systemd-run", "--version"], text=True).splitlines()[0],
               "probe_results": records, "scientific_resource_gate_passed": False,
               "scope": "synthetic instrumentation qualification only"}
    (outdir / "probe.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("Verified descendant CPU, simultaneous memory peak, fresh counter reset, failed-command and orphan rejection")


def validate_probes(records):
    serial, parallel = records["serial"], records["parallel"]
    if any(not r["measurement_complete"] or r["descendants_remaining"]
           or r["cpu_seconds"] < .4 for r in (serial, parallel)):
        raise ValueError("CPU accounting missed descendants or measurement incomplete")
    # Two simultaneous allocations must exceed sequential reuse. Broad bounds
    # allow interpreter/kernel overhead without accepting cumulative child peaks.
    if not (30e6 < serial["peak_memory_bytes"] < 75e6
            and 65e6 < parallel["peak_memory_bytes"] < 140e6
            and parallel["peak_memory_bytes"] > serial["peak_memory_bytes"] + 20e6):
        raise ValueError("memory counters do not distinguish simultaneous from serial allocations")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--worker", choices=("serial", "parallel", "failure", "orphan"), help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(args.outdir, args.worker)
    else:
        run_probe(args.outdir)
