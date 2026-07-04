# Higher-Order Convection Discretization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Replace Ember's 1st-order upwind convection discretization with a
config-selectable, limited, (mostly) 2nd-order upwind scheme plus a
trapezoidal continuity march, validated by example baselining and
grid-convergence studies.

**Architecture:** A new shared derivative kernel (`ConvectionDifferencer`)
used by both `ConvectionSystemUTW` and `ConvectionSystemY`; config plumbing
follows the existing `splittingMethod` pattern; legacy scheme retained
bit-identical behind `firstOrderUpwind`. See the approved spec:
`docs/superpowers/specs/2026-07-04-convection-discretization-design.md` —
**the spec is the authority for all numerical formulas (§3) and architecture
(§4). Read it before starting any task.**

**Tech Stack:** C++17 (Eigen via `dvec`), gtest (`test/*.cpp`), Python/Cython
config bridge (`python/ember/input.py`, `_ember.pxd/.pyx`), SCons via pixi.

## Plan style note (intentional deviation from skill defaults)

Per project owner instruction, this plan is **deliberately high-level**: tasks
specify files, interfaces, requirements, and verification — not line-level
code. Implementing agents are expected to make code-level decisions
themselves, guided by the spec. Where the spec pins down a formula or name,
it is binding; everything else is implementer's choice following existing
codebase patterns.

## Global constraints

- Config option: `General.convectionScheme`, values `"secondOrderLimited"`
  (default) and `"firstOrderUpwind"` (legacy). Kernel enum:
  `ConvectionDifferencer::Scheme::{FirstOrderUpwind, SecondOrderLimited}`.
- `firstOrderUpwind` must reproduce current behavior **bit-identically**
  (including the forced forward difference at `j=0` and the rectangle-rule
  continuity march).
- Limiter: van Albada (spec §3.1). Boundary slopes `σ[0] = σ[jj] = 0`.
  Trapezoidal continuity march per spec §3.3, gated on the scheme option.
- Each `ConvectionSystemUTW`/`ConvectionSystemY` instance owns its own
  differencer (TBB thread safety).
- Quasi-2D: fix the uninitialized-`v` upwind bug by passing node-wise `vr`
  velocities to the kernel (spec §3.4). Behavior change accepted; quasi-2D is
  an orphan configuration and NOT a merge blocker.
- Untouched subsystems: `DiffusionSystem`, `TridiagonalIntegrator`,
  `SourceSystem`, `SplitSolver`, `updateCrossTerms`, grid adaptation.
- Build: new `src/*.cpp` and `test/*.cpp` files are auto-globbed by SConstruct
  — no build-system edits needed. Build: `pixi run build`. Tests:
  `pixi run test` (runs gtest suite and Python tests).
- Commit after each task with message prefix `convection:` and the task ID,
  e.g. `convection: [1.1] add ConvectionDifferencer kernel`.

## Resumption protocol (for interrupted sessions / fresh contexts)

1. Read the spec, then this plan. Checkboxes below track completed steps;
   `git log --oneline --grep=convection:` shows completed task commits.
2. Baselines and convergence artifacts live in `test/convergence/` (created in
   Task 0.1); its README documents regeneration.
3. If a task is half-done (commits exist but boxes unchecked or vice versa),
   trust the git history + test results over the checkboxes, fix the
   checkboxes, and continue.
4. Delegation notes per task suggest the least-capable adequate model
   (Sonnet for mechanical/pattern-following work, Opus for numerics-sensitive
   work). The orchestrating agent reviews all subagent output between tasks.

---

## Phase 0 — Baselines (no code changes)

### Task 0.1: Example baseline harness and capture

**Delegation:** Sonnet (scripting; no numerics judgment needed).

**Files:**
- Create: `test/convergence/README.md` (how to regenerate everything here)
- Create: `test/convergence/run_baselines.py`
- Create: `test/convergence/baselines/*.json` (captured outputs)

**Interfaces:**
- Produces: `run_baselines.py --outdir test/convergence/baselines [--cases ...]`
  which runs the curated example subset and writes one JSON per case
  containing: case name, config summary (including scheme option once it
  exists), final time, grid size N, key scalars (peak T, consumption/flame
  speed where defined, integral heat release), and the final `x`, `T`, and
  major-species profiles.
- Produces: `test/convergence/compare_baselines.py baselineA baselineB`
  reporting relative differences in the scalars and profile L2/L∞ norms
  (interpolated to a common grid), with a pass/fail threshold argument.

**Requirements:**
- Curated cases (from spec §6.3): `example_single`, `example_diffusion`,
  `example_twin`, `example_cylindrical_outward`, `example_cylindrical_inward`,
  `example_laminarFlameSpeed`.
- Runtimes may be shortened (reduced `tEnd`, smaller mechanisms are already
  used) as long as each case runs long enough to be a meaningful indicator;
  record every deviation from the stock example config in the JSON.
- Baselines MUST be captured from the unmodified pre-implementation
  code (current `main`, commit recorded in the JSON). If Phase 1 has already
  landed when this runs, capture from a clean checkout of the recorded
  pre-implementation commit instead.

**Steps:**
- [x] Write `run_baselines.py` and `compare_baselines.py` with README
- [x] Run the six cases; verify each completes and writes plausible JSON
      (no NaNs, T within physical bounds)
- [x] Re-run one case twice to quantify run-to-run reproducibility (thread
      scheduling noise); record the observed tolerance floor in the README —
      this becomes the comparison threshold for Task 1.5
- [x] Commit artifacts: `convection: [0.1] example baseline harness + baselines`

---

## Phase 1 — Convection upgrade

### Task 1.1: `ConvectionDifferencer` kernel + unit tests

**Delegation:** Opus (numerics-sensitive; formulas fully specified in spec §3.1–3.2).

**Files:**
- Create: `src/convectionDifferencer.h`, `src/convectionDifferencer.cpp`
- Test: `test/test_convectionDifferencer.cpp` (gtest, follows existing
  `test/test_diffusionSystem.cpp` conventions)

**Interfaces (binding, from spec §4.1):**
```cpp
class ConvectionDifferencer {
public:
    enum class Scheme { FirstOrderUpwind, SecondOrderLimited };
    void setScheme(Scheme s);
    void resize(size_t nPoints);   // sizes σ scratch
    //! dy/dx at nodes 0..jj-1; upwind branch per-node from sign(v[j]).
    //! Right-boundary node jj is handled by callers.
    void computeDerivatives(const dvec& y, const dvec& v,
                            const OneDimGrid& grid, dvec& dydx);
};
```
- `FirstOrderUpwind`: bit-identical to the legacy loops in
  `ConvectionSystemUTW::f()` / `ConvectionSystemY::f()` (forward difference
  when `v[j] < 0 || j == 0`, else backward).
- `SecondOrderLimited`: van Albada σ pass + upwind face-reconstruction pass
  per spec §3.1; `σ[0] = σ[jj] = 0`; `dydx[0] = 0` in this path; asserts
  `nPoints >= 3`.

**Test requirements (TDD — write tests first, watch them fail):**
- Exactness: constant field → identically 0; linear field → exact derivative
  (both schemes, mixed-sign velocity).
- Bit-identical legacy check: `FirstOrderUpwind` vs. a copy of the legacy
  loop, random data and velocities, exact equality.
- Order verification: smooth profile (e.g. tanh) on randomized nonuniform
  grids with adjacent-spacing ratios up to a `uniformityTol`-like factor
  (~2); observed convergence order of `SecondOrderLimited` ≥ 2 in interior
  L2 norm under grid refinement. (Build small standalone `OneDimGrid`
  instances; see how `test_diffusionSystem.cpp` constructs grids and call
  `setSize`/`updateValues`.)
- No-new-extrema: for monotone data, reconstructed face values stay within
  the data range (test via the derivative's sign structure or by exposing
  face values internally — implementer's choice, but the property must be
  tested).
- Upwind-branch correctness: velocity field with sign changes → derivative at
  each node uses the correct side (compare against hand-computed small case).

**Steps:**
- [x] Write failing unit tests covering all requirements above
- [x] Implement kernel; iterate until tests pass (`pixi run test`)
- [x] Commit: `convection: [1.1] ConvectionDifferencer kernel + unit tests`

### Task 1.2: Config plumbing for `convectionScheme`

**Delegation:** Sonnet (mechanical; mirrors `splittingMethod` exactly).

**Files:**
- Modify: `python/ember/input.py` (General options, near `splittingMethod`)
- Modify: `python/ember/_ember.pxd`, `python/ember/_ember.pyx`
- Modify: `src/readConfig.h`, `src/readConfig.cpp`

**Interfaces:**
- Produces: `ConfigOptions::convectionScheme` (`std::string`), values
  `"secondOrderLimited"` (default) / `"firstOrderUpwind"`.
- Python: `General.convectionScheme = StringOption("secondOrderLimited",
  ("firstOrderUpwind",), level=2)` with docstring stating: 2nd-order in
  smooth regions, local 1st-order limiter fallback at extrema;
  `firstOrderUpwind` is the pre-1.7 legacy scheme.

**Requirements:**
- Follow the `splittingMethod` pattern end-to-end (option definition →
  pxd/pyx bridge → ConfigOptions member). Validation of allowed values
  happens via `StringOption`'s existing machinery; C++-side parsing to the
  kernel enum happens in Task 1.3 and must throw on unknown strings.

**Steps:**
- [x] Add option across the four plumbing layers
- [x] Rebuild (`pixi run build`); verify `Config(General(convectionScheme=
      'firstOrderUpwind')).validate()` passes and an invalid value fails
      validation (quick Python check, can be a throwaway snippet)
- [x] Commit: `convection: [1.2] convectionScheme config option`
      (note: `readConfig.cpp` needed no change — string options bridge via
      `.pyx` directly, matching `splittingMethod`)

### Task 1.3: Integrate kernel into `ConvectionSystemUTW` + trapezoidal march

**Delegation:** Opus (touches continuity logic and stagnation-point branches).

**Files:**
- Modify: `src/convectionSystem.h`, `src/convectionSystem.cpp`
  (`ConvectionSystemUTW::f()`, `ConvectionSystemSplit::setTolerances()`,
  `ConvectionSystemSplit::resize()`)

**Interfaces:**
- Consumes: `ConvectionDifferencer` (Task 1.1), `ConfigOptions::convectionScheme`
  (Task 1.2).
- Produces: `ConvectionSystemSplit` parses the option string to
  `ConvectionDifferencer::Scheme` once in `setTolerances()` (throw
  `DebugException` on unknown value) and propagates it to `utwSystem` and,
  in Task 1.4, to species systems.

**Requirements:**
- Replace the inline `dTdx/dUdx/dWdx` upwind loop with three kernel calls
  using velocity `rV` for branch selection (unchanged sign convention).
- Trapezoidal `rV` march per spec §3.3, in all continuity-BC branches
  (forward and backward marches), gated on `SecondOrderLimited`;
  `FirstOrderUpwind` keeps the rectangle rule verbatim.
- Boundary rows (`j=0`, `j=jj`) and all `dUdt/dTdt/dWdt` assembly unchanged.
- The differencer is a member of `ConvectionSystemUTW`, resized alongside the
  system.

**Verification:**
- [ ] Full test suite passes (`pixi run test`)
- [ ] With `firstOrderUpwind`: `test/python/test_flame_configs.py` results
      unchanged from pre-change behavior (spot-check one config by diffing
      final T profile against a pre-change run, or rely on Task 1.5 baseline
      comparison if preferred — state which was done in the commit message)
- [ ] Commit: `convection: [1.3] UTW system on ConvectionDifferencer + trapezoidal continuity`

### Task 1.4: Integrate kernel into `ConvectionSystemY` (incl. quasi-2D fix)

**Delegation:** Opus (quasi-2D velocity semantics need care).

**Files:**
- Modify: `src/convectionSystem.h`, `src/convectionSystem.cpp`
  (`ConvectionSystemY::f()`, `ConvectionSystemSplit::configureSolver()`)

**Interfaces:**
- Consumes: `ConvectionDifferencer`, scheme enum from `ConvectionSystemSplit`
  (Task 1.3).

**Requirements:**
- Replace the inline per-node upwind logic with one kernel call per `f()`
  evaluation; each `ConvectionSystemY` owns its differencer (TBB safety).
- Standard path: velocity = interpolated `v` (as today).
- Quasi-2D path: build a node-wise advecting-velocity array from
  `vrInterp->get(x[j], t)` and pass it to the kernel — this fixes the
  uninitialized-`v` upwind-direction bug (spec §3.4). The `ydot` assembly
  keeps its existing quasi-2D form (`-vr * dYdx / vz`).
- Boundary rows (`j=0` inflow balance, `j=jj` outflow) unchanged.

**Verification:**
- [ ] Full test suite passes with both scheme settings
- [ ] Commit: `convection: [1.4] species convection on ConvectionDifferencer + quasi-2D velocity fix`

### Task 1.5: Regression + baseline comparison

**Delegation:** Sonnet (run, compare, report; escalate anomalies).

**Files:**
- Create: `test/convergence/baselines-phase1/*.json` (post-change runs)
- Modify: `test/convergence/README.md` (record results summary)

**Requirements:**
- Run the Task 0.1 case set with `firstOrderUpwind`: must match Phase 0
  baselines within the reproducibility tolerance recorded in Task 0.1.
- Run with `secondOrderLimited`: expect small, physically consistent shifts
  (report the deltas; no pass/fail threshold, but flag anything anomalous —
  sign changes in flame speed trends, unbounded Y, NaNs, CVODE step-count
  blowup vs. legacy — for review by the orchestrator/user before merge).
- Record CVODE step counts (`getNumSteps()` totals are logged when
  `debugTimesteps` is on) for both schemes on at least `example_single` and
  `example_twin`.

**Steps:**
- [ ] Run both scheme settings across the case set; write comparison report
      into README (tables of scalar deltas + step counts)
- [ ] Verify `firstOrderUpwind` parity; investigate/fix any mismatch before
      proceeding (bit-identical expectation modulo thread noise)
- [ ] Commit: `convection: [1.5] phase-1 regression + baseline comparison`

**Phase 1 exit criteria:** all unit + regression tests pass with both
schemes; `firstOrderUpwind` parity confirmed; `secondOrderLimited` results
reviewed with no anomalies. Phase 1 is independently mergeable here.

---

## Phase 2 — Convergence studies and diffusion decision gate

### Task 2.1: Convergence study scripts

**Delegation:** Sonnet (scripting to a fixed protocol), Opus review of the
protocol implementation.

**Files:**
- Create: `test/convergence/run_convergence.py`
- Modify: `test/convergence/README.md`

**Interfaces:**
- Produces: `run_convergence.py --case {strained,twin,cylindrical} --scheme
  {firstOrderUpwind,secondOrderLimited} [--damp-const X]` writing per-run
  JSON (N, grid tolerances, consumption speed, peak T, CVODE step totals,
  wall time) into `test/convergence/results/`.
- Produces: `plot_convergence.py` generating error-vs-N plots (reference =
  finest `secondOrderLimited` run per case) for both schemes.

**Requirements (protocol from spec §6.4):**
- Case A: unbounded strained premixed H₂ flame (base on `example_single` /
  test configs; FixedValue + ZeroGradient BCs).
- Case B: twin flame (ControlVolume BC) and cylindrical flame (α=1). Both
  BC paths must appear in the study; if runtime forces a choice, twin is
  required, cylindrical strongly preferred.
- Resolution ladder via `vtol`/`dvtol`/`gridMax` sweeps (≥5 rungs spanning
  ~4× in N); fixed strain/composition per case; steady-state termination.

**Steps:**
- [ ] Write scripts; smoke-test one rung per case
- [ ] Commit: `convection: [2.1] convergence study scripts`

### Task 2.2: Run studies, analyze, and decide Phase 3

**Delegation:** Runs by Sonnet; analysis and written findings by Opus;
**decision requires user review — do not proceed into Phase 3 without it.**

**Files:**
- Create: `test/convergence/results/` artifacts +
  `docs/superpowers/specs/2026-07-04-convection-discretization-design.md`
  gets a short "Phase 2 findings" addendum section (append, don't rewrite).

**Requirements:**
- Produce error-vs-N curves for both schemes, all cases; report observed
  orders and the grid-size reduction at matched error.
- `dampConst` relaxation trial on one case with `secondOrderLimited` to
  quantify additional coarsening headroom (spec §6.4).
- CVODE step-count comparison (limiter smoothness check, spec §5).
- Written recommendation on the Phase 3 gate: does diffusion error dominate
  at target resolutions? Evidence: if `secondOrderLimited` error-vs-N
  flattens toward 2nd-order-diffusion-limited behavior while convection
  refinement no longer helps, the gate opens.

**Steps:**
- [ ] Run full matrix; generate plots and findings addendum
- [ ] Present findings + Phase 3 recommendation to user; record decision in
      the addendum
- [ ] Commit: `convection: [2.2] convergence study results + phase-3 decision`

---

## Phase 3 — Diffusion upgrade (CONDITIONAL — only if Phase 2 gate opens)

Deliberately not planned in detail (spec §7). If the gate opens:

- Option (a) minimal: harmonic-mean face diffusivity within the existing
  3-point tridiagonal structure (`DiffusionSystem::get_A`).
- Option (b) full: 5-point 4th-order conservative flux stencils +
  pentadiagonal BDF integrator (extension of `TridiagonalIntegrator`).

**Process requirement:** run a fresh brainstorming/design pass (new spec
section or standalone spec) scoped by the Phase 2 evidence, then extend this
plan with concrete tasks. Do not improvise Phase 3 directly from this
paragraph.

---

## Future work (out of scope, documented)

Conservative flux-form framework (Approach C): see spec Appendix A. The
kernel's face-reconstruction internals are the intended reuse seam
(`ỹ_{j±1/2}` values); keep them cleanly separable when implementing Task 1.1,
but do not add speculative API surface for it (YAGNI).
