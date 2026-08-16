# repeat-binder-gonogo — control run and results

Validates a binder-scoring pipeline (Boltz-2, via `proto-tools` on Modal)
against repeat-protein complexes with known outcomes, per
`bench/repeat-binder-gonogo/task.md`. This is a control run — nothing here
designs a binder; it decides whether the scoring stack can be trusted to.

## Read this first: which result is "the" result

There are two runs in here, and only one of them is a submission.

| | `bench/repeat-binder-gonogo/results/` | `demo/full-sweep/` |
|---|---|---|
| Command | `python -m gonogo --scorer boltz2 run` | `python -m gonogo run --ignore-veto` |
| Follows the task's cost policy | yes | **no, by explicit override** |
| G2 (specificity) | **fail** | fail (same numbers) |
| G1 / G3 / G4 | `not_run` (task rule: G2 failed → stop) | G1/G3 still `not_run`, **G4 pass** |
| Verdict | this is the real answer | exploratory — "what would G1/G3/G4 have said anyway" |

**The finding: G2 fails.** DARPin 4m3 scores 78.7 against mouse cathepsin B
(the true 65.7 nM binder) and 39.7 against human cathepsin B (the true
non-binder, 83% identical) — a 38.9-point margin in the right direction. But
Boltz-2's own seed-to-seed noise at 3 seeds is *wider* than that margin
(52.7-point spread on B4 alone: one seed scored it 67, another 14.6). Per the
task's rule — margin must exceed the model's own uncertainty band — that is a
fail, not a low pass. **The pipeline is not yet certified for specificity
claims**, and per the task's cost policy, nothing past G2 should be bought
until that's fixed.

`demo/full-sweep/` exists because we wanted to know what G1/G3/G4 would have
said anyway, as a diagnostic. Its own trajectory stamps this explicitly:
*"veto overridden ... This is a demo sweep, not a submission."* Both runs fail
the bench's pytest verifier identically, because both correctly report G2 as
failed — `bench/repeat-binder-gonogo/verifier/test.sh` scores a failed or
`not_run` G2 as reward 0 by design; that is the gate the whole validation
exists for.

## Two follow-up diagnostics, not part of either run

- **`bench/repeat-binder-gonogo/results/g1_dockq.json`** — G1 is written in
  DockQ, which isn't in the `proto-tools` catalogue. This installs the
  standalone `DockQ` package directly and runs it against all 8 Tier C
  structures. Result: every comparison fails inside DockQ's own PDB parser
  before producing a score. So even with DockQ available, G1 doesn't close for
  free — something about how the predicted structures are written trips
  DockQ's loader, and that's the next thing to debug, not a config toggle.
- **`bench/repeat-binder-gonogo/results/g2_seed_depth.json`** — resamples B1
  and B4 out to 12 seeds and asks whether more seeds would flip G2. Under the
  harness's own band (max − min across seeds), the answer is no — it fails all
  the way to n=12. Under a Welch's-t confidence interval, it flips to pass at
  n≥4. Marked `"does not amend the shipped verdict"` on purpose: this is
  evidence about which **band definition** is right, not permission to
  overturn G2 after the fact.

## Layout

```
gonogo/                       the harness (see gonogo/README.md)
bench/repeat-binder-gonogo/   the task: cases.yaml, REFERENCES.md, the verifier,
                               and results/ (the real run + the two diagnostics above)
demo/selftest/                harness self-test, synthetic scorer, meaningless numbers —
                               exercises the gate arithmetic and verifier wiring only
demo/full-sweep/               the --ignore-veto exploratory run described above
docs/gonogo-feasibility.md    the full write-up: what's feasible against task.md,
                               what's blocked and why, cost breakdown, five things
                               to fix in cases.yaml
scripts/proto_tools_smoke_test.py   smoke-tests a deployed proto-tools GPU service
```

## Rerunning

```bash
python -m gonogo preflight          # free: what can run, what it costs
python -m gonogo inputs             # free: resolve cases, re-derive cases.yaml's own claims
python -m gonogo --scorer boltz2 run
python -m gonogo verify
```

Needs `boltz2-prediction` deployed on your own Modal account
(`proto-tools deploy --apps boltz2 --env proto-env`) — billed to whoever runs
it. `python -m gonogo --scorer synthetic --out demo/selftest run --allow-synthetic`
exercises the harness itself without a GPU or a bill.

## Next steps, in order

1. Fix `gonogo`'s scorer to run more seeds per complex (5+), or switch the G2
   band to a CI rather than a range — `g2_seed_depth.json` is the evidence for
   deciding which.
2. Debug the DockQ/Boltz-2 structure mismatch in `g1_dockq.json`; once G1 can
   run, it should run over Tier A and Tier C per the task.
3. Only after G2 passes for real: buy Tier A (63 variants, ~$13) and finish
   G3/G4 over the full case set.
