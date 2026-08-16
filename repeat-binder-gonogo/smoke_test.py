"""Smoke-test the three GPU_NO_HIGH_END services against Modal.

    .venv/bin/python backups/proto-tools-gpu-patch-2026-08-15/smoke_test.py

One representative tool per patched app, run with the tool's own canonical
example input so nothing here depends on local data. Reports where each call
actually ran and how long it took — the first call per app pays a container
start and a weight download, so a slow first result is expected.
"""

from __future__ import annotations

import sys
import time

from proto_tools.mcp.tools import run_tool

ALL_CASES = [
    ("esm2-embedding", "esm2"),
    ("esmfold-prediction", "esmfold"),
    ("proteinmpnn-sample", "proteinmpnn"),
    ("boltz2-prediction", "boltz2"),
]

# Optional positional args select a subset by app slug, e.g. `... smoke_test.py boltz2`.
wanted = set(sys.argv[1:])
CASES = [case for case in ALL_CASES if not wanted or case[1] in wanted or case[0] in wanted]
if not CASES:
    raise SystemExit(f"no case matches {sorted(wanted)}; known apps: {sorted({app for _, app in ALL_CASES})}")

failures = 0
for tool_key, app in CASES:
    print(f"\n=== {tool_key} (app: {app}) ===", flush=True)
    started = time.time()
    try:
        result = run_tool(tool_key, use_example=True, device="modal")
    except Exception as exc:  # a failed smoke test is a result, not a crash
        failures += 1
        print(f"  FAILED after {time.time() - started:.1f}s: {type(exc).__name__}: {exc}")
        continue
    elapsed = time.time() - started
    if not result.get("ok", True):
        failures += 1
        print(f"  FAILED after {elapsed:.1f}s: {result.get('error')}")
        continue
    payload = result.get("result", {})
    summary = {
        key: (f"<{len(value)} items>" if isinstance(value, (list, dict)) else value)
        for key, value in payload.items()
        if key in {"tool_id", "execution_time", "success", "warnings", "errors"}
    }
    print(f"  ok in {elapsed:.1f}s  ran_on={result.get('ran_on')}  {summary}")
    if result.get("saved_files"):
        print(f"  files: {result['saved_files']}")

print(f"\n{len(CASES) - failures}/{len(CASES)} smoke tests passed")
raise SystemExit(1 if failures else 0)
