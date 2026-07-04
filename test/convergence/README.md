# Convergence / regression baselines

This directory captures reference outputs from a curated subset of Ember's
example configurations, to be diffed against runs of the modified convection
solver later in this project (see the convection-scheme design doc and
implementation plan, Phase 1 onward).

## Contents

- `run_baselines.py` — runs the curated cases and writes one JSON file per
  case into `baselines/`.
- `compare_baselines.py` — compares two baseline JSON files (scalars and
  profile L2/L-infinity norms on a common interpolated grid) against a
  pass/fail threshold.
- `baselines/*.json` — the captured baseline outputs (committed to git).

## Curated cases

Six cases from spec §6.3, each a faithful transcription of the corresponding
`python/ember/examples/example_*.py` script (only `Paths.outputDir`/`logFile`
are redirected into a scratch directory — no physics changes):

- `example_single`
- `example_diffusion`
- `example_twin`
- `example_cylindrical_outward`
- `example_cylindrical_inward`
- `example_laminarFlameSpeed`

No case needed its `tEnd`/tolerance shortened relative to the stock example:
each one's own early-termination criterion (steady heat-release or dT/dt
threshold) kicks in well before its stock `tEnd` cap, so wall-clock time per
case is already small (8-50 s; ~2.3 minutes for all six). Any future
deviation from a stock example configuration must be recorded both in the
case-builder function in `run_baselines.py` and in that case's
`deviations_from_stock` JSON field.

## Provenance

Baselines in this repository were captured from commit
`5552b7d203291c1e62f4f95ca9627ae26b95d853` (tip of `convection-scheme` at the
time Task 0.1 ran), which at that point contained only planning docs — no
solver code changes — so it stands in for "unmodified pre-implementation
code." Each JSON's `commit` field records the exact commit this rule applies
to (with a `-dirty` suffix if tracked files differed from that commit;
untracked scratch files are ignored). If Phase 1 has already landed on this
branch by the time you need a fresh Phase-0-equivalent baseline, regenerate
from a clean checkout of that same recorded commit instead of `HEAD`.

## Regenerating

```
pixi run python test/convergence/run_baselines.py \
    --outdir test/convergence/baselines \
    [--cases example_single example_diffusion ...] \
    [--workdir build/test/baselines-work] \
    [--retries 3]
```

Run from the repository root. `--workdir` is scratch space for Ember's own
HDF5/log output (gitignored under `build/`); only the summarized JSON in
`--outdir` is meant to be committed.

Note: with the stock `nThreads` values (2 or 4 depending on case), the
solver occasionally (rarely) aborts mid-run with `CVODE Integrator had too
many errors` — this was observed once in roughly a handful of attempts on
`example_single`, and is apparently a manifestation of the thread-scheduling
nondeterminism referenced in the task brief, not something introduced by
this harness. `run_baselines.py` retries a failed case up to `--retries`
times (default 3) before giving up; each JSON records how many attempts it
took in its `attempts` field.

## Comparing two baselines

```
pixi run python test/convergence/compare_baselines.py \
    test/convergence/baselines/example_single.json \
    /path/to/other/example_single.json \
    --threshold 0.02
```

Prints relative differences for each scalar and normalized L2/L-infinity
norms for each profile (T and each captured major species, linearly
interpolated onto a shared grid spanning the overlap of the two x-domains).
Exits 0 if every relative difference/norm is at or below `--threshold`, 1
otherwise — usable as a scripted pass/fail gate (e.g. for Task 1.5).

## Reproducibility floor (thread-scheduling noise)

Per the task brief, one case (`example_single`, the case that uses the
stock `nThreads=4` and happened to exhibit the CVODE failure above) was run
three times from the same commit/config and compared pairwise with
`compare_baselines.py`:

| Comparison             | worst relative scalar diff | worst profile L2 | worst profile L-inf |
|-------------------------|----------------------------|-------------------|----------------------|
| run1 vs run2            | 5.2e-7 (heat release)       | 1.4e-7 (Y_OH)     | 6.7e-7 (Y_CO)        |
| run1 vs run3            | 2.5e-7 (heat release)       | 2.3e-7 (Y_OH)     | 1.4e-6 (Y_OH)        |
| run2 vs run3            | 7.7e-7 (heat release)       | 1.6e-7 (Y_O)      | 1.3e-6 (Y_O)         |

All three runs used identical configuration and commit; the ~1e-6-level
differences are consistent with ordinary floating-point summation-order
noise from multi-threaded execution, not a meaningful physical difference.

**Observed tolerance floor: ~1e-6 relative (scalars and profile norms).**
Any Phase 1 comparison should treat differences at or below roughly 1e-5 as
run-to-run noise, not evidence of a behavior change. The `compare_baselines.py`
default threshold (`0.02`, i.e. 2%) is set well above this floor, with large
margin to also tolerate small, legitimate discretization differences from a
new convection scheme while still catching real regressions. **The other,
more significant reproducibility caveat is qualitative, not quantitative:**
under the stock multi-threaded configuration, a run can occasionally fail
outright (CVODE integrator error) rather than produce a slightly different
answer. Task 1.5 tooling should be prepared to retry a failed run rather
than treat a crash as a silent baseline mismatch.

## JSON schema

Each `baselines/<case>.json` contains:

- `case` — example name
- `commit` — git commit SHA the run was captured from (see Provenance above)
- `generated_at_utc` — capture timestamp
- `config_summary` — `Config.stringify()` text of the full effective configuration
- `deviations_from_stock` — list of human-readable notes on anything that
  differs from running the stock example script directly
- `scheme` — convection discretization scheme option; `null` at this commit
  because that option does not exist yet (Phase 0). Once Phase 1 adds a
  `scheme` option to `General`, this field should be populated.
- `attempts` — number of attempts `run_baselines.py` needed for this case (see above)
- `runtime_seconds` — wall-clock time for the successful attempt
- `final_time` — solver's `tNow` at termination
- `grid_size` — number of grid points at termination
- `scalars` — `peak_T`, `consumption_speed`, `heat_release_rate_integral`, `flame_position`
  (any scalar that isn't physically meaningful for that case's geometry is `null`;
  see `scalar_notes` for why)
- `scalar_notes` — explanations for any `null`/discarded scalar (e.g.
  `consumption_speed` is ill-defined whenever the domain's two boundary
  temperatures are nearly equal, which happens for `example_single`, a
  premixed flame opposed by a cold inert at the same temperature as the
  reactants, and for `example_diffusion`, whose fuel and oxidizer streams
  are both preheated to the same temperature)
- `profiles` — final `x`, `T`, and mass-fraction profiles for up to 8 "major"
  species (peak mass fraction >= 1e-3 over the profile, most-abundant first)
