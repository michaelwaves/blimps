#!/bin/bash
# Verifier: run the pytest suite, turn pass-rate into a reward in [0, 1].
# Exit 0 after writing reward.txt; nonzero means verifier infrastructure failure.
set -uo pipefail

REPORT=/logs/verifier/report.json
mkdir -p /logs/verifier

python -m pytest "$(dirname "$0")/test_outputs.py" \
  -p no:randomly -q \
  --json-report --json-report-file="$REPORT" \
  >/logs/verifier/pytest.log 2>&1

if [ ! -f "$REPORT" ]; then
  echo "pytest produced no report; see /logs/verifier/pytest.log" >&2
  echo "0.0" > /logs/verifier/reward.txt
  exit 0
fi

python - "$REPORT" <<'PY'
import json, sys

summary = json.load(open(sys.argv[1]))["summary"]
passed = summary.get("passed", 0)
total = summary.get("total", 0) - summary.get("skipped", 0)
reward = passed / total if total else 0.0

# G2 is the gate the whole validation exists for: failing it caps the run.
for test in json.load(open(sys.argv[1]))["tests"]:
    if "g2_separates" in test["nodeid"] and test["outcome"] == "failed":
        reward = 0.0
        break

open("/logs/verifier/reward.txt", "w").write(f"{reward:.4f}\n")
print(f"reward={reward:.4f} ({passed}/{total} checks)")
PY
