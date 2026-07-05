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
- `total_convection_steps` — sum, over every global timestep in the run, of
  `ConvectionSystemSplit::getNumSteps()` (see Task 1.5 section below for how
  this is derived). `null` if `Debug.timesteps` was off or no log file was
  found.

## Phase 1 regression results (Task 1.5)

This section records the results of running the curated case set with each
of the two `General.convectionScheme` values and comparing against the
Phase 0 baselines above. Raw JSONs live in `baselines-phase1/*.json`
(filenames suffixed `_firstOrderUpwind` / `_secondOrderLimited`; a couple of
extra `_firstOrderUpwind_run2`/`_run3` repeatability-check files for
`example_single` are also included, see the caveat at the end of this
section).

### Harness extension: `--scheme`

`run_baselines.py` gained a `--scheme {firstOrderUpwind,secondOrderLimited}`
flag (default: no override, i.e. whatever each case's config defaults to,
currently `secondOrderLimited`). When given, it sets
`conf.general.convectionScheme.value` for every case in the run before
`conf.validate()`/`conf.evaluate()`, and the resulting `deviations_from_stock`
list and `scheme` JSON field record the override. This required no changes to
`compare_baselines.py`.

`run_baselines.py` also now sums the per-global-timestep convection CVODE
step counts that `FlameSolver::finishStep()` writes to the log file whenever
`Debug.timesteps` is enabled (the default: `Debug.timesteps = True`), e.g.
lines of the form:

```
t = 0.000100 (dt = 1.000e-04) [C: 12]
```

The trailing integer is `ConvectionSystemSplit::getNumSteps()`, which is
**reset every global timestep** (`FlameSolver::setState()` reinitializes the
split convection CVODE integrators once per step), so it is a per-step count,
not a running total. `total_convection_steps()` in `run_baselines.py` sums
every occurrence across the run's log file to produce a whole-run total,
recorded in the new `total_convection_steps` JSON field (kept outside the
`scalars` dict so it is not swept into `compare_baselines.py`'s pass/fail
scalar check).

### Commands used

```
pixi run python test/convergence/run_baselines.py \
    --outdir test/convergence/baselines-phase1 \
    --workdir build/test/baselines-work-phase1 \
    --scheme firstOrderUpwind --suffix _firstOrderUpwind

pixi run python test/convergence/run_baselines.py \
    --outdir test/convergence/baselines-phase1 \
    --workdir build/test/baselines-work-phase1-solim \
    --scheme secondOrderLimited --suffix _secondOrderLimited
```

All twelve runs (six cases x two schemes) completed on the first attempt
(no `CVODE Integrator had too many errors` retries needed for the main
sweep; one retry was needed for an incidental third repeatability run of
`example_single`, consistent with the previously-documented flake rate).

### `firstOrderUpwind` parity vs. Phase 0

Compared with `compare_baselines.py --threshold 0.02` (the tool's
pass/fail gate) against `test/convergence/baselines/*.json`:

| Case                          | Worst relative diff (scalar or profile norm) | Where            | Result |
|-------------------------------|-----------------------------------------------|-------------------|--------|
| example_single                | 1.50e-4                                        | Y_O L-inf         | PASS   |
| example_diffusion              | 1.30e-9                                        | flame_position    | PASS   |
| example_twin                  | 1.31e-5                                        | Y_O2 profile      | PASS   |
| example_cylindrical_outward    | 1.92e-7                                        | Y_O2 profile      | PASS   |
| example_cylindrical_inward     | 7.84e-8                                        | Y_OH profile      | PASS   |
| example_laminarFlameSpeed      | 2.06e-7                                        | Y_OH profile      | PASS   |

All six cases pass with wide margin under the 0.02 gate; five of six are at
or below the ~1e-6 "reproducibility floor" recorded in the Phase 0 README.
`example_single` is the exception — see the caveat below.

**Caveat — `example_single` termination-time sensitivity (not scheme-specific):**
While investigating why `example_single`'s worst diff (1.5e-4) was ~100x the
1e-6 floor (still comfortably inside the 0.02 gate, but conspicuous), two
extra `firstOrderUpwind` repeats of the *identical* config/commit were run
(`example_single_firstOrderUpwind_run2.json`, `_run3.json`). Their
`final_time` at termination differed noticeably from the first run and from
each other (0.0144 / 0.0186 / 0.0166 s), and pairwise profile comparisons
between these same-code repeats showed differences up to ~14% (Y_OH
L-infinity), i.e. *larger* than the scheme-comparison numbers above. This
case's default termination criterion (`TerminationCondition(measurement='Q')`,
a tight RMS-steadiness check on heat release rate) is evidently sensitive
enough to thread-scheduling floating-point noise that it can trigger at
noticeably different simulated times from run to run — consistent with the
occasional `CVODE Integrator had too many errors` abort already documented
for this exact case in the Phase 0 README (one of the three extra runs hit
this on its first attempt and succeeded on retry). This is a pre-existing
property of the case's termination logic, reproduced identically under
`firstOrderUpwind` alone (no scheme comparison involved), so it is not
evidence of a defect introduced by the convection-scheme change — but it does
mean the ~1e-6 floor documented in Phase 0 does not hold reliably for this
specific case on this machine/build, and any single-run comparison of
`example_single` (in either direction) carries substantially more inherent
noise than the other five cases. Flagged here for visibility; not blocking
per the 0.02 pass/fail gate, but worth keeping in mind for Phase 2
convergence-study work that reuses this case.

### `secondOrderLimited` deltas vs. Phase 0 (no pass/fail threshold)

Scalar relative differences:

| Case                       | peak_T   | consumption_speed | heat_release_integral | flame_position |
|-----------------------------|----------|--------------------|-------------------------|-----------------|
| example_single               | 1.08e-3  | n/a (null)         | 3.34e-2                 | 7.39e-2         |
| example_diffusion             | 1.41e-3  | n/a (null)         | 4.95e-3                 | 2.06e-2         |
| example_twin                 | 1.72e-3  | 8.00e-3            | 6.15e-3                 | 1.76e-2         |
| example_cylindrical_outward   | 4.53e-3  | 3.09e-3            | 8.91e-3                 | 1.87e-2         |
| example_cylindrical_inward    | 2.43e-3  | 1.68e-2            | 1.73e-2                 | 1.50e-3         |
| example_laminarFlameSpeed     | 1.06e-3  | 2.98e-2            | 3.20e-2                 | 1.22e-3         |

Profile L2 / L-infinity norms (worst species shown; full output in the
`compare_baselines.py` invocations below):

| Case                       | T (L2 / Linf)        | Worst species (L2 / Linf)      |
|-----------------------------|------------------------|----------------------------------|
| example_single               | 0.034 / 0.099          | Y_CO2: 0.039 / 0.116             |
| example_diffusion             | 0.0036 / 0.0099        | Y_OH: 0.0030 / 0.0159            |
| example_twin                 | 0.028 / 0.119          | Y_OH: 0.052 / 0.339              |
| example_cylindrical_outward   | 0.014 / 0.051          | Y_OH: 0.031 / 0.146              |
| example_cylindrical_inward    | 0.0043 / 0.018         | Y_H2: 0.0095 / 0.037             |
| example_laminarFlameSpeed     | 0.011 / 0.084          | Y_O: 0.023 / 0.251               |

All deltas are consistent with a small, physically-reasonable shift from
adding second-order limited reconstruction: peak temperature moves by
<0.5%, heat release rate integral and consumption speed move by 1-3%
(no sign changes — see below), and the largest pointwise (L-infinity) norms
are concentrated in thin radical layers (OH, O), which is expected: a small
absolute shift in flame position/thickness translates into a large
*normalized, pointwise* error for a species whose profile is a sharp,
narrow spike, even when the underlying physical change is modest (the L2
norms, which average over the whole domain rather than taking a single
worst point, are 3-10x smaller than the L-infinity values in every case).
`example_twin` also has a different major-species set at the 1e-3 peak-Y
cutoff (`O` in Phase 0 vs. `NO` under `secondOrderLimited`), reflecting a
small shift in trace-radical peak levels rather than a missing/spurious
species.

**Anomaly checks (per task brief):**
- Sign changes in flame-speed trend: none. `consumption_speed` stays
  positive and within ~3% of Phase 0 for all four cases where it's
  physically meaningful (`example_twin`, `example_cylindrical_outward`,
  `example_cylindrical_inward`, `example_laminarFlameSpeed`).
- Unbounded/negative mass fractions: none found (checked `min(Y) >= 0` for
  every captured species profile in every `secondOrderLimited` run).
- NaNs / non-finite values: none (`run_baselines.py` raises on any
  non-finite `T`/`Y`; no run raised, and `min(T) >= 300 K` / sane `max(T)`
  in all six cases).
- CVODE step-count blowup: none. See step-count table below — the largest
  increase is 1.35x (`example_single`), well short of a "blowup."

### CVODE convection step counts

Whole-run totals (`total_convection_steps`, see harness extension above),
for all six cases (task requires at least `example_single`/`example_twin`,
included here for all six for completeness):

| Case                       | firstOrderUpwind | secondOrderLimited | ratio (2nd/1st) |
|-----------------------------|-------------------|----------------------|------------------|
| example_single               | 803,865           | 1,082,716            | 1.35             |
| example_diffusion             | 162,713           | 164,043              | 1.01             |
| example_twin                 | 142,045           | 150,427              | 1.06             |
| example_cylindrical_outward   | 143,580           | 168,666              | 1.18             |
| example_cylindrical_inward    | 227,462           | 283,964              | 1.25             |
| example_laminarFlameSpeed     | 431,872           | 559,069              | 1.30             |

`secondOrderLimited` consistently takes modestly more convection CVODE steps
(1.0-1.35x) than `firstOrderUpwind`, plausibly reflecting the added
stiffness/nonlinearity of the van Albada flux limiter. No case shows an
order-of-magnitude increase.

### Reproducing these comparisons

```
pixi run python test/convergence/compare_baselines.py \
    test/convergence/baselines/example_single.json \
    test/convergence/baselines-phase1/example_single_firstOrderUpwind.json \
    --threshold 0.02

pixi run python test/convergence/compare_baselines.py \
    test/convergence/baselines/example_single.json \
    test/convergence/baselines-phase1/example_single_secondOrderLimited.json \
    --threshold 1.0   # no meaningful gate for this comparison; raises the
                       # threshold purely to suppress "EXCEEDS THRESHOLD" noise
```

### Continuity-BC smoke tests (final whole-branch review fix)

All six curated cases above (and the Task 1.5 sweep) use the default
`General.continuityBC='fixedLeft'`. The `stagnationPoint` and
`fixedTemperature` continuity boundary conditions
(`ContinuityBoundaryCondition::Zero` / `::Temp` in `convectionSystem.cpp`)
take a structurally different path through `ConvectionSystemUTW::f()` — the
`rV` march runs forward *and* backward from an interior `jContBC` rather than
forward from `rV[0]` — and, per the final whole-branch review, had never
actually executed under the new `secondOrderLimited` trapezoidal march before
this check (nor, for that matter, under `firstOrderUpwind`, since no test or
baseline configures anything other than `fixedLeft`).

This was a does-it-execute-sanely smoke check, not a physics/convergence
study: four short runs (`continuityBC` in `{stagnationPoint,
fixedTemperature}` x `convectionScheme` in `{secondOrderLimited,
firstOrderUpwind}`), each based on `TestPremixedStrained` from
`test/python/test_flame_configs.py` (H2/O2/Ar premixed counterflow flame,
`h2o2.yaml` mechanism, `nThreads=1`) with `TerminationCondition(tEnd=0.006,
measurement=None)` — at the default `globalTimestep=2e-5 s` this is ~300
global timesteps, enough for the flame to develop a real interior stagnation
point and temperature gradient (both BCs need one to lock onto) while
staying fast. `fixedTemperature` requires `General.splittingMethod
='balanced'`, which is already Ember's default, so no extra config was
needed. Script: `test/convergence/smoke_continuity_bc.py` (`pixi run python
test/convergence/smoke_continuity_bc.py`).

Pass criteria: run completes without exception, finite T/Y throughout, and
peak T within 250-3000 K.

| continuityBC       | convectionScheme    | Completed | final_time | grid_size | peak_T (K) | min_T (K) | finite |
|--------------------|---------------------|-----------|------------|-----------|------------|-----------|--------|
| stagnationPoint     | secondOrderLimited   | yes       | 0.006      | 65        | 1592.5     | 300.0     | yes    |
| stagnationPoint     | firstOrderUpwind     | yes       | 0.006      | 60        | 1516.7     | 300.0     | yes    |
| fixedTemperature    | secondOrderLimited   | yes       | 0.006      | 68        | 1557.2     | 300.0     | yes    |
| fixedTemperature    | firstOrderUpwind     | yes       | 0.006      | 67        | 1538.2     | 300.0     | yes    |

All four runs completed on the first attempt (no CVODE retries needed; the
whole matrix was also re-run once in full as a repeatability check with
identical qualitative results). No NaN/Inf in any T or Y profile; all peak
temperatures are well inside the 250-3000 K physical bound. Qualitative
scheme comparison: `secondOrderLimited` vs. `firstOrderUpwind` peak-T
differs by ~5.0% for `stagnationPoint` and ~1.2% for `fixedTemperature` —
larger than the ~0.1-0.5% seen in the fully-converged Task 1.5 sweep above,
but expected here: these are short, still-transient runs (fixed `tEnd`
rather than a steady-state termination criterion), so a modest shift in
front position/speed from the different convection discretization is
amplified relative to a converged profile. Neither BC nor either scheme
shows a wild divergence, a sign flip, or any instability symptom (grid size,
`min_T`, and `finite` are all consistent across the matrix), so both
previously-unexercised branches are confirmed to execute sanely under both
convection schemes. No numerics bug found; this closes the residual test-gap
item from the final whole-branch review.
