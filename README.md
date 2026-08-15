# blimps

Harbor tasks for high-risk-bio evals.

## Setup

```bash
uv sync
source .venv/bin/activate   # or prefix commands with `uv run`
harbor --version
proto-tools doctor
```

`harbor[modal]` and `proto-tools` are declared in `pyproject.toml` and locked in
`uv.lock` (proto-tools pinned to an exact commit of
`github.com/michaelwaves/proto-tools`), so every clone resolves the same
versions. Credentials still come from your own environment (`modal token new`,
`harbor auth login`, `ANTHROPIC_API_KEY`).

proto-tools caches model weights and per-tool environments under `PROTO_HOME`
(default `~/.proto`); set that variable to move them off your home volume.

## Running a task

```bash
T=tasks/ppi-design

uv run harbor task start-env -p $T -e modal          # interactive env
uv run harbor run -p $T -a oracle -e modal           # reference solution
uv run harbor run -p $T -a claude-code -m claude-opus-5 -e modal
```
