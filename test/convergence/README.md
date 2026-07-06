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
- `smoke_continuity_bc.py` — does-it-execute-sanely smoke check for the
  `stagnationPoint`/`fixedTemperature` continuity BCs (see "Continuity-BC
  smoke tests" section below).
- `run_convergence.py` / `plot_convergence.py` — the Task 2.2 grid-convergence
  study harness (Task 2.1 builds these scripts and smoke-tests them; Task 2.2
  runs the full matrix). See "Grid-convergence study (Task 2.1/2.2)" below.
- `results/` — output directory for `run_convergence.py` JSON files. Empty
  (placeholder `.gitkeep` only) as of Task 2.1; the full study's results are
  committed here in Task 2.2.

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

## Grid-convergence study (Task 2.1/2.2)

This section documents `run_convergence.py` and `plot_convergence.py`, the
harness for the grid-convergence study required by spec §6.4. **Task 2.1
(this task) writes and smoke-tests these scripts; it does not run the full
study matrix.** Task 2.2 runs the complete (case x scheme x rung) sweep and
commits the results.

### Cases

Three cases, chosen to cover both convection-BC families named in spec §6.4
(`ConvectionSystemUTW::f()`'s two branches over `grid.leftBC`) with a fixed
strain rate and composition per case (only the grid tolerances vary across
rungs):

| Case | `--case` | Composition / strain | Left BC | Right BC | Geometry |
|---|---|---|---|---|---|
| A | `strained` | H2:1.0 / O2:1.0,AR:4.0, phi=0.3, a=800/s | `ZeroGradient` | `FixedValue` | planar, alpha=0 |
| B1 | `twin` | CH4:1.0/air, phi=0.70, a=100/s | `ControlVolume` | `FixedValue` | planar/twin, alpha=0 |
| B2 | `cylindrical` | CH4:0.5,H2:0.5/air, phi=0.60, a=500/s | `ControlVolume` | `FixedValue` | cylindrical, alpha=1 |

- **`strained`** (Case A, "unbounded"): `General(unburnedLeft=False,
  fixedBurnedVal=False)` on the default planar geometry. Tracing
  `FlameSolver::updateBC()` (src/flameSolver.cpp:619-638) and
  `OneDimGrid::updateBoundaryIndices()` (src/grid.cpp:839-846): with
  `unburnedLeft=False`, `jb=0` (burned/product index) and `ju=jj`
  (unburned/reactant index); with `fixedBurnedVal=False` and this being
  neither a twin nor a cylindrical flame, the left-boundary branch falls
  through every special case (not `WallFlux`, not `ju==0`, not `jb==0 &&
  fixedBurnedVal`, not twin/cylindrical centering) to the `else`:
  `BoundaryCondition::ZeroGradient` — an open, unclamped product-side
  boundary that the adaptive grid can extend outward as needed ("unbounded").
  The right boundary (`ju==jj`, the physical reactant inlet) stays
  `BoundaryCondition::FixedValue`. This is the literal `FixedValue`/
  `ZeroGradient` pairing spec §6.4 names for Case A, and it is structurally
  distinct from Case B's `ControlVolume` inflow-balance path. Composition and
  strain rate follow `TestPremixedStrained`
  (test/python/test_flame_configs.py) and `smoke_continuity_bc.py`.
- **`twin`** (Case B1): `General(twinFlame=True, unburnedLeft=False)`, which
  resolves to `BoundaryCondition::ControlVolume` at the symmetry/stagnation
  plane (x=0). Physical setup follows `example_twin.py` /
  `case_twin` in `run_baselines.py`.
- **`cylindrical`** (Case B2): `General(flameGeometry='cylindrical',
  unburnedLeft=False, fixedLeftLocation=True)`, i.e. curvature parameter
  `alpha=1` in `OneDimGrid` (vs. `alpha=0` for the planar cases above). Also
  resolves to `ControlVolume`, exercising the curved-geometry (`rphalf`
  weighting) path instead of the planar-twin path. Physical setup follows
  `example_cylindrical_outward.py` / `case_cylindrical_outward` in
  `run_baselines.py`.

All three use `General(nThreads=1)` for run-to-run determinism (the
baseline harness's thread-scheduling noise, documented above, is an
unwanted confound for a convergence-order study). `strained` uses the small
`h2o2.yaml` mechanism (fast); `twin`/`cylindrical` use the stock `gri30.yaml`
mechanism matching their source examples. Termination is measurement-based
in all three cases (steady state, not a fixed `tEnd`): `strained` and
`cylindrical` use `TerminationCondition(measurement='dTdt')` explicitly;
`twin` relies on the default `measurement='Q'` (as in `example_twin.py`,
which only overrides `tEnd` as a hard cap). `tEnd` is set generously in each
case as a safety cap, not the expected termination trigger.

### Resolution ladder

Six rungs (`RUNGS` in `run_convergence.py`), index 0 (coarsest) to 5
(finest), sweeping `Grid(vtol, dvtol, gridMax)` together (each rung also
scales `gridMax` down so the max-spacing cap doesn't become the binding
constraint before `vtol`/`dvtol` do):

| rung | vtol | dvtol | gridMax |
|---|---|---|---|
| 0 | 0.24 | 0.40 | 4.0e-4 |
| 1 | 0.16 | 0.27 | 2.5e-4 |
| 2 | 0.11 | 0.18 | 1.6e-4 |
| 3 | 0.075 | 0.12 | 1.0e-4 |
| 4 | 0.050 | 0.080 | 6.3e-5 |
| 5 | 0.033 | 0.055 | 4.0e-5 |

`vtol`/`dvtol` shrink by a factor of ~1.5x per rung (~7.3x coarsest to
finest). This is a **nominal** ladder: run_convergence.py defines >= 5 rungs
intended to span roughly a 4x range in N per spec §6.4, but the actual N
achieved by each rung is configuration-dependent (a thin, highly-strained
flame needs more points per unit `vtol` than a slow one) and was not
verified end-to-end for all rungs under Task 2.1 (only rung 0 was
smoke-tested per case, see below). **Task 2.2 should confirm the achieved N
values actually span roughly 4x-8x and retune the `RUNGS` table (values
only, keep 6 rungs) if a case's ladder is off**, e.g. if two adjacent rungs
happen to produce the same N (redundant) or the span comes out much
narrower/wider than intended.

`--damp-const X` overrides `Grid.dampConst` (default 7) for a single run,
independent of the rung ladder — this is spec §6.4's "trial with relaxed
`dampConst` to quantify the achievable grid coarsening", run as a one-off
comparison against the matching (case, scheme, rung) run without the
override, not swept across all rungs.

### Usage

```
# List the rung ladder
pixi run python test/convergence/run_convergence.py --list-rungs

# Run every rung for one (case, scheme) pair (Task 2.2 usage)
pixi run python test/convergence/run_convergence.py \
    --case strained --scheme secondOrderLimited

# Run a single rung (e.g. for a smoke test)
pixi run python test/convergence/run_convergence.py \
    --case strained --scheme secondOrderLimited --rung 0

# Run an explicit subset of rungs
pixi run python test/convergence/run_convergence.py \
    --case twin --scheme firstOrderUpwind --rungs 0 2 4

# One-off relaxed-dampConst trial
pixi run python test/convergence/run_convergence.py \
    --case cylindrical --scheme secondOrderLimited --rung 3 --damp-const 15
```

Full option list: `--case {strained,twin,cylindrical}` and `--scheme
{firstOrderUpwind,secondOrderLimited}` (both required unless `--list-rungs`
is given), `--damp-const X`, `--rung N` / `--rungs N [N ...]` (default: all
6 rungs), `--outdir` (default `test/convergence/results`), `--workdir`
(default `build/test/convergence-work`), `--retries` (default 3, same
retry-on-exception behavior as `run_baselines.py` — see the known-flake note
below).

Output: one JSON file per run,
`<outdir>/<case>_<scheme>_rung<N>[_damp<X>].json`.

### JSON schema (`run_convergence.py` output)

- `case`, `scheme`, `rung`, `damp_const` (`null` unless `--damp-const` given)
- `commit`, `generated_at_utc`, `config_summary` — provenance, same
  convention as `run_baselines.py`
- `grid_tolerances` — `vtol`, `dvtol`, `gridMax`, `gridMin`, `dampConst`, read
  back from the concrete (post-`evaluate()`) config, i.e. the tolerances
  actually used, not just the rung table's nominal values
- `N` — final grid size (`len(solver.x)`)
- `scalars.consumption_speed`, `scalars.peak_T` — the two convergence
  metrics required by spec §6.4; `consumption_speed` is `null` (with a note
  in `scalar_notes`) if the domain's two boundary temperatures are nearly
  equal, mirroring `run_baselines.py`'s guard (not expected for any of these
  three cases' physical setups, but checked defensively)
- `total_convection_steps` — whole-run sum of the per-global-timestep
  `ConvectionSystemSplit::getNumSteps()` log values, via the same
  `total_convection_steps()`/regex helper `run_convergence.py` imports from
  `run_baselines.py`
- `runtime_seconds` — wall-clock time for the run
- `final_time` — solver's `tNow` at termination
- `attempts` — number of attempts needed (see retry note below)

### Plotting

```
pixi run python test/convergence/plot_convergence.py \
    [--cases strained twin cylindrical] \
    [--resultsdir test/convergence/results] \
    [--outdir test/convergence/results/plots]
```

For each case and metric (`consumption_speed`, `peak_T`), loads every
non-`--damp-const` run, takes the **largest-N `secondOrderLimited` run as
the reference**, and plots log-log relative error vs. N for both schemes on
one figure (`<case>_<metric>.png`). Skips (with a message, not an error) any
case/metric that doesn't have at least a reference run plus 2 comparable
points for a scheme, so it can be run against a partial results directory
mid-sweep.

### Known flake (inherited from run_baselines.py)

Same as documented above for the baseline harness: the solver occasionally
aborts with `CVODE Integrator had too many errors` under stiff/multi-step
runs. `run_convergence.py` retries a failed run up to `--retries` times
(default 3, matching `run_baselines.py`) before giving up; each JSON records
`attempts`.

### Task 2.1 smoke test

One rung (rung 0, the coarsest) per case, default scheme
(`secondOrderLimited`), run via:

```
pixi run python test/convergence/run_convergence.py \
    --case strained --scheme secondOrderLimited --rung 0 \
    --outdir build/test/convergence-smoke --workdir build/test/convergence-work-smoke
pixi run python test/convergence/run_convergence.py \
    --case twin --scheme secondOrderLimited --rung 0 \
    --outdir build/test/convergence-smoke --workdir build/test/convergence-work-smoke
pixi run python test/convergence/run_convergence.py \
    --case cylindrical --scheme secondOrderLimited --rung 0 \
    --outdir build/test/convergence-smoke --workdir build/test/convergence-work-smoke
```

(Output redirected to `build/` — already gitignored — rather than
`test/convergence/results/`, since a single coarsest-rung run per case isn't
a representative "result" of the study; `results/` stays an empty
placeholder for Task 2.2's full sweep.)

| case | scheme | rung | completed | N | final_time | peak_T (K) | consumption_speed (m/s) | convection_steps | wall time (s) |
|---|---|---|---|---|---|---|---|---|---|
| strained | secondOrderLimited | 0 | yes | 55 | 0.0046 | 1558.8 | 0.3569 | 16,245 | 0.9 |
| twin | secondOrderLimited | 0 | yes | 76 | 0.0074 | 1842.6 | 0.1702 | 130,626 | 13.2 |
| cylindrical | secondOrderLimited | 0 | yes | 64 | 0.0082 | 1790.9 | 0.2260 | 144,926 | 13.7 |

All three runs completed on the first attempt (no retries needed). Every
JSON has all required fields populated (`N`, `grid_tolerances`, both
scalars, `total_convection_steps`, `runtime_seconds`); no NaN/Inf values
anywhere in any of the three JSON files (checked recursively over every
float field, not just the headline scalars); `peak_T` is physically
sane (1550-1850 K, consistent with the lean/moderate-equivalence-ratio
mixtures used); `consumption_speed` is positive and O(0.1-0.4 m/s), as
expected for these mixtures/strain rates. This confirms the harness
produces valid output for all three cases and both BC families
(`ZeroGradient`/`FixedValue` for `strained`, `ControlVolume`/`FixedValue`
for `twin` and `cylindrical`) at the coarsest rung; the full multi-rung,
multi-scheme sweep (needed to actually assess convergence order) is Task
2.2's job.

## errTol-based grid adaptation (Tasks G1-G7)

Tasks G1-G6 replaced the `vtol`/`dvtol` grid-refinement criterion with a
dimensional, scheme-aware local-error budget (`Grid.errTol`, default
`1e-4`; `vtol`/`dvtol` are now deprecated/ignored, see
`python/ember/input.py`). Design and full evidence chain:
`docs/superpowers/specs/2026-07-05-error-based-grid-adaptation-design.md`
(addendum §A covers the G6 errTol-ladder acceptance run in detail).
Calibration of the `secondOrderLimited` error coefficient
(`errCoeff`, `1/15 -> 0.0139`) and the chosen default `errTol`:
`test/convergence/results/calibration-notes.md`. Raw ladder JSONs/plots
(untracked, regenerate via `run_convergence.py`/`plot_convergence.py`):
`test/convergence/results-errtol-final/`. The tol -> measured-QoI-error
table referenced throughout the design doc is
`test/convergence/results-errtol-final/plots/errtol_error_table.md`
(machine-generated by `plot_convergence.py`).

### Grid-settling analyzer (`analyze_settling.py`)

`analyze_settling.py <workdir> [<workdir> ...]` inspects the numbered
`profNNNNNN.h5` outputs in a `run_convergence.py`/`run_baselines.py`
work directory and classifies the grid-size trajectory as `SETTLED`,
`NOT SETTLED`, or `TOO FEW OUTPUTS` (fewer than 4 profile snapshots,
not judgeable). Used to distinguish a healthy overshoot-then-plateau
grid history from the pre-errTol P2.4 failure signature (monotonic,
never-settling grid growth).

`grid_settled(traj, tail_frac=0.25, tol_pts=None)` calls a trajectory
settled if the grid size varies by at most `tol_pts` over the last
`tail_frac` of outputs. `tol_pts=None` (the default) uses a *relative*
tolerance, `max(3, ceil(0.05 * median(tail)))`, rather than a fixed
point count: an earlier fixed `tol_pts=3` produced false `NOT SETTLED`
verdicts on 14/36 runs of the G6 errTol ladder once N grew into the
hundreds (an absolute point tolerance doesn't scale with N). Verified
against the ladder data: `strained`/`secondOrderLimited` rung 5
(`build/test/convergence-work-errtol-final/strained_secondOrderLimited_rung5`,
tail spread 7 pts at tail-median N=218) now reports `SETTLED`, while the
old pre-errTol P2.4 failure case
(`build/test/convergence-work/strained_secondOrderLimited_rung5`, tail
spread 94 pts, monotonic growth to N=369) still correctly reports
`NOT SETTLED`. (A first attempt at a 1%-of-median coefficient, floated
as an example in the design addendum §A.5, was checked against this
same data and found to still round down to the old absolute floor of 3
for the settled case — 5% was adopted instead, with margin above the
~3.2% the verified case needs and well below the ~32% relative spread
of the genuine failure case.) Pass an explicit numeric `tol_pts` to
restore the old fixed-tolerance behavior.

### Task G7: baseline regression at errTol defaults

With `errTol=1e-4` now the library default (previously the six curated
cases below all ran under `vtol=0.12`/`dvtol=0.2`), Task G7 reran the
Phase-1 curated case set at library defaults (no `--scheme` override,
i.e. whatever each stock example config already specifies — five of six
cases default to `secondOrderLimited`) and compared against the
Phase-1 `_secondOrderLimited` baselines. Outputs:
`test/convergence/results/baselines-errtol/` (this run's JSONs).

```
pixi run -- python test/convergence/run_baselines.py \
    --outdir test/convergence/results/baselines-errtol \
    --workdir build/test/baselines-work-errtol

pixi run -- python test/convergence/compare_baselines.py \
    test/convergence/baselines-phase1/<case>_secondOrderLimited.json \
    test/convergence/results/baselines-errtol/<case>.json \
    --threshold 1.0   # compare_baselines' pass/fail gate isn't
                       # calibrated for this cross-default comparison;
                       # read the printed deltas directly instead
```

**Result: 5/6 cases regress cleanly; `example_single` fails to complete
at the new default.**

| Case | peak_T rel. diff | consumption_speed rel. diff | N (phase-1 -> errTol default) |
|---|---|---|---|
| example_cylindrical_inward | 6.07e-5 | 1.24e-3 | 146 -> 170 |
| example_cylindrical_outward | 1.91e-5 | 7.50e-4 | 125 -> 154 |
| example_diffusion | 2.50e-5 | n/a (null, near-equal boundary temperatures) | 185 -> 212 |
| example_laminarFlameSpeed | 5.68e-4 | 5.56e-4 | 158 -> 183 |
| example_twin | 5.00e-6 | 5.62e-5 | 154 -> 169 |
| example_single | **did not complete** (`CVODE Integrator had too many errors`, deterministic, see below) | | 194 -> (n/a) |

All five completing cases pass the amended spec §4 scalar criterion with
wide margin (worst delta 5.68e-4, i.e. 0.057%, well under the 0.5%
threshold). **However, `N` grew for every one of the five** (+10% to
+23%), which the spec §4 "same or fewer" clause does not license — see
the ledger close-out entry for the open question this raises about the
`errTol=1e-4` default (already pending owner sign-off from G5/G6 for
other reasons).

**`example_single` regression:** the stock config (`nThreads=4`,
`disc` geometry, cold-inert counterflow) raised
`CVODE Integrator had too many errors` on all 6 attempts (`--retries 6`),
including with `nThreads=1` (rules out the previously-documented
thread-scheduling flake — this is deterministic). Bisecting `errTol`
directly on this case (`nThreads=1`, otherwise stock config): `errTol=5e-4`
and `errTol=2e-4` both complete cleanly (N=152, N=197); `errTol=1e-4`
(the library default) fails every time. This case was not part of the
G1-G6 three-case (`strained`/`twin`/`cylindrical`) ladder study, so this
robustness gap at the exact default value was not previously exercised.
Per instruction, the default was **not** retuned to paper over this — it
is recorded here as a finding for owner review alongside the other
`errTol=1e-4` default questions.
