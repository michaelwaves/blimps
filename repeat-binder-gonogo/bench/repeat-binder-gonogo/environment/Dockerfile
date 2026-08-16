FROM ubuntu:24.04

RUN apt-get update -qq \
 && apt-get install -y -qq curl python3 python3-pip python3-venv \
 && rm -rf /var/lib/apt/lists/*

# Verifier dependencies. The agent's own scientific stack (proto-tools et al.)
# is supplied by the runner, not baked in here.
RUN python3 -m pip install --break-system-packages --no-cache-dir \
      pytest pytest-json-report pyyaml

WORKDIR /app

RUN mkdir -p /logs/verifier /logs/agent /logs/artifacts /app/results
