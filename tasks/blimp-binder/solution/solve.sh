#!/bin/bash
set -euo pipefail

mkdir -p /logs/outputs
/app/foundry/.venv/bin/python /solution/solve.py
