# Design: Higher-order convection discretization for Ember

**Date:** 2026-07-04
**Status:** Approved design, pending implementation
**Scope:** Convection operator (advective derivatives + continuity march), with a
conditional later phase for the diffusion operator. Conservative flux-form
rewrite (Approach C) documented as future work in the Appendix.

## 1. Context and motivation

Ember solves the 1D flame equations with Strang/balanced operator splitting
(`SplitSolver`), integrating convection, diffusion, and chemical source terms
separately. The spatial discretization is a legacy of fully-coupled 1D solvers,
where keeping the Jacobian tightly banded forces compact low-order stencils.
The splitting removes that constraint: the convection sub-problem is integrated
explicitly per-component (CVODE/Adams), so wider or nonlinear (limited)
stencils carry almost no structural cost.

Current discretization:

- **Convection** (`convectionSystem.cpp`): 1st-order upwind for all advective
  derivatives (`U`, `T`, `Wmx`, each species `Y_k`); continuity (`rV`) marched
  with a 1st-order rectangle rule from the anchor point `jContBC`.
- **Diffusion** (`diffusionSystem.cpp`): conservative 3-point flux-divergence
  scheme, 2nd-order, tridiagonal, integrated with a custom BDF1/BDF2 solver
  (Thomas algorithm).
- **Cross terms** (Soret, correction flux, enthalpy flux): 2nd-order
  face-based, explicit, frozen over each split step.

The 1st-order upwind convection is the dominant spatial error source. The grid
module's `dampVal` mechanism caps grid spacing based on the
convective/diffusive ratio specifically to suppress the numerical diffusion
this scheme generates — i.e., the current scheme forces grid refinement to
hide its own error.

## 2. Goals and non-goals

**Primary goal:** accuracy per grid point — obtain the same solution quality
with substantially fewer grid points by eliminating the O(h) numerical
diffusion of 1st-order upwind convection.

**Goals:**

- Replace 1st-order upwind advective derivatives with a limited, upwind-biased,
  (mostly) 2nd-order scheme valid on the strongly nonuniform adaptive grid.
- Upgrade the continuity (`rV`) march from rectangle to trapezoidal quadrature.
- Make the scheme selectable via the config system; new scheme is the default;
  legacy scheme retained bit-identical as `firstOrderUpwind`.
- Validate with grid-convergence studies covering all boundary-condition
  paths: unbounded strained flame **and** twin/cylindrical flames
  (ControlVolume BC, α=1 geometry, stagnation-point continuity anchor).

**Non-goals (this phase):**

- Diffusion operator changes — deferred to a conditional Phase 3, decided by
  the Phase 2 convergence experiments.
- Conservative flux-form reformulation (Approach C) — see Appendix.
- Grid adaptation changes (`dampVal`/`dampConst` retuning is a Phase 2
  experiment, not a code change).
- Restoring/validating the quasi-2D path (known orphan; see §7).

## 3. Numerical scheme

### 3.1 Limited upwind advective derivative (interior nodes `j = 1..jj−1`)

One-sided slopes at node `j`:

```
s⁻ = (y[j] − y[j−1]) / hh[j−1]
s⁺ = (y[j+1] − y[j]) / hh[j]
```

Limited node slope, **van Albada** limiter (default):

```
σ[j] = (s⁻ s⁺)(s⁻ + s⁺) / ((s⁻)² + (s⁺)²)   if s⁻·s⁺ > 0
σ[j] = 0                                      otherwise
```

van Albada is chosen over minmod/MC because it is smooth away from slope sign
changes, which matters because CVODE/Adams evaluates the RHS many times per
step and its step-size/order selection benefits from RHS smoothness.
(MC-limited weighted-central slope is a documented alternative if sharper
front resolution is ever needed.)

Advective derivative by upwind-biased face reconstruction:

```
∂y/∂x |_j ≈ (ỹ_{j+1/2} − ỹ_{j−1/2}) / dlj[j]
```

with reconstructions from the upwind side of each face:

- For `v[j] > 0`:  `ỹ_{j+1/2} = y[j]   + σ[j]·hh[j]/2`,
  `ỹ_{j−1/2} = y[j−1] + σ[j−1]·hh[j−1]/2`
- For `v[j] < 0`:  `ỹ_{j+1/2} = y[j+1] − σ[j+1]·hh[j]/2`,
  `ỹ_{j−1/2} = y[j]   − σ[j]·hh[j−1]/2`

Properties (verify in unit tests):

- 2nd-order at the node **regardless of grid nonuniformity**, given 2nd-order
  slopes σ (Taylor expansion: the h² terms cancel in the face difference).
- Constant fields → derivative exactly 0. Linear fields → exact derivative
  (limiter returns the exact slope when `s⁻ = s⁺`).
- Falls back to 1st-order upwind exactly where the limiter clips (extrema),
  which is what boundedness of `Y ∈ [0,1]` requires.
- Stencil: 4 points (`j−2..j+1` for `v>0`, mirrored for `v<0`).

### 3.2 Closures and special points

- `σ[0] = σ[jj] = 0`: reconstructions at the outermost faces fall back to the
  adjacent cell value (locally 1st-order, monotone). These are far-field
  regions with small gradients.
- Boundary nodes `j=0` and `j=jj` are **unchanged**: ControlVolume/WallFlux
  inflow balance, FixedValue/ZeroGradient suppression, and the 1st-order
  one-sided outflow at `jj` all stay as-is.
- Upwind branch selection per node: `sign(rV[j])` in `ConvectionSystemUTW`,
  `sign(v[j])` in `ConvectionSystemY`, exactly as today. The `j==0` forced
  forward difference in the legacy loop is preserved in the legacy path.
- Stagnation point (`Zero`/`Temp` continuity BCs): no special handling beyond
  the per-node branch; `|v| → 0` there so the advective term vanishes. The
  `rV` linear seeding around `xVzero` is unchanged.
- Minimum grid size: scheme is well-defined for `nPoints ≥ 3`; kernel asserts.

### 3.3 Trapezoidal continuity march

Replace (current form, rectangle rule with face radius and left-node integrand):

```
rV[j+1] = rV[j] − hh[j] · rphalf[j] · (drhodt[j] + rho[j]·β·U[j])
```

with the trapezoidal rule using node radii:

```
S[j]    = drhodt[j] + rho[j]·β·U[j]
rV[j+1] = rV[j] − hh[j]/2 · (r[j]·S[j] + r[j+1]·S[j+1])
```

applied consistently to the forward and backward marches from `jContBC` in all
continuity-BC branches. At the axis (α=1, `r[0]=0`) the integrand vanishes
naturally.

The trapezoidal march is **tied to the same scheme option** as the derivative
upgrade: `firstOrderUpwind` keeps the rectangle rule so the legacy path stays
bit-identical. (Accepted trade-off: A/B comparisons bundle the derivative and
march upgrades.)

### 3.4 Quasi-2D bug fix (behavior change)

In `ConvectionSystemY::f()` with `quasi2d == true`, `update_v()` is never
called, so the upwind test `if (v[j] < 0)` reads **uninitialized memory** and
the upwind direction is arbitrary. The new kernel takes the advecting velocity
explicitly; the quasi-2D branch fills a node-wise velocity array from
`vrInterp` and passes it in, fixing the bug as a side effect. Quasi-2D results
may change. Per project decision, quasi-2D is an orphaned configuration:
behavior changes there are acceptable and it is not a merge blocker (§7).

## 4. Architecture

**Approach chosen:** shared derivative kernel ("Approach B") — one small
component used by both convection system classes. Rejected alternatives:
in-place duplicated stencils in each `f()` ("Approach A", violates DRY across
4–5 call sites); conservative flux-form framework ("Approach C", deferred —
see Appendix).

### 4.1 New component: `ConvectionDifferencer`

New files `src/convectionDifferencer.h` / `.cpp`:

```cpp
class ConvectionDifferencer
{
public:
    enum class Scheme { FirstOrderUpwind, SecondOrderLimited };

    void setScheme(Scheme s);
    void resize(size_t nPoints);   // sizes the σ scratch vector

    //! Advective derivative dy/dx at nodes 0..jj-1, upwind branch chosen
    //! per-node from sign(v[j]). Right-boundary node handled by callers.
    void computeDerivatives(const dvec& y, const dvec& v,
                            const OneDimGrid& grid, dvec& dydx);
};
```

- `FirstOrderUpwind` path reproduces the existing loops **bit-identically**,
  including the forced forward difference at `j=0`.
- `SecondOrderLimited` path: one pass to compute σ, one pass for the face
  differences. `dydx[0]` is set to 0 in this path (unused by callers).
- Unit-testable in isolation against analytic functions on synthetic
  nonuniform grids.

### 4.2 Call sites

- `ConvectionSystemUTW::f()`: three kernel calls (`U`, `T`, `Wmx`) with
  velocity `rV`, replacing the inline upwind loop. Trapezoidal `rV` march
  implemented in-place in `f()` (continuity logic, not a derivative),
  gated on the scheme.
- `ConvectionSystemY::f()`: one kernel call with `v` (quasi-2D: node-wise
  `vr` array from `vrInterp`). `ydot` assembly, BC rows, and `splitConst`
  handling untouched.
- **Threading:** each `ConvectionSystemUTW` / `ConvectionSystemY` instance
  owns its differencer (species systems run in parallel under TBB; per-system
  scratch avoids sharing). Cost: two scratch `dvec`s per system.

### 4.3 Config plumbing (follows `splittingMethod` pattern)

- `python/ember/input.py`:
  `General.convectionScheme = StringOption("secondOrderLimited", ("firstOrderUpwind",), level=2)`
  with a docstring noting the limited scheme is 2nd-order in smooth regions
  with local 1st-order fallback at extrema.
- `python/ember/_ember.pxd` / `_ember.pyx`: declare and assign the string.
- `src/readConfig.h`: `std::string convectionScheme;` member.
- `ConvectionSystemSplit::setTolerances(options)` (already receives
  `ConfigOptions`, already called from `FlameSolver::initialize`) parses the
  string to the enum once; `resize()` / `configureSolver()` propagate it to
  `utwSystem` and each `speciesSystems[k]`. Unknown strings throw at startup.

### 4.4 Explicitly untouched

`DiffusionSystem`, `TridiagonalIntegrator`, `SourceSystem`, `SplitSolver`,
cross-term computation (`updateCrossTerms`), and all grid adaptation logic.

## 5. Edge cases and risks

| Case | Handling |
|---|---|
| Extrema / steep fronts | Limiter clips σ → local 1st-order fallback; `Y` boundedness preserved; `correctMassFractions()` remains as backstop |
| Stagnation point | Per-node upwind branch as today; term vanishes as `\|v\|→0`; `xVzero` seeding unchanged |
| `j=1` needs `σ[0]`; `j=jj−1` needs `σ[jj]` | Defined as 0 — no out-of-bounds stencil |
| Small grids (`nPoints` ≈ 4–5) | Well-defined for `nPoints ≥ 3`; kernel asserts |
| CVODE/Adams step-size behavior vs. limiter nonsmoothness | van Albada is C¹ away from sign changes; monitor `getNumSteps()` in A/B runs — step-count blowup is a pre-merge red flag |
| Regression drift | Existing Python tests anchor to physical values (Cantera equilibrium); expected to pass under new default; verify both schemes |
| Quasi-2D behavior change | Expected (bug fix §3.4); quasi-2D is an orphan configuration with no automated tests — acceptable, not a merge blocker |

## 6. Testing and validation

### 6.1 Unit tests (`test/test_convectionDifferencer.cpp`, gtest)

- Exactness on constant and linear fields (both schemes).
- Observed 2nd-order convergence of `SecondOrderLimited` on smooth profiles
  over randomized nonuniform grids (spacing ratios up to `uniformityTol`).
- No-new-extrema property on monotone data.
- Bit-identical match of `FirstOrderUpwind` against a copy of the legacy loop.
- Upwind-branch correctness for mixed-sign velocity fields.

### 6.2 Regression tests

- Existing C++ (gtest) and Python (`test/python/`) suites pass with **both**
  scheme settings.

### 6.3 Example-based baselines (pre-implementation)

The unit test suite is known to be incomplete; the examples in
`python/ember/examples/` serve as additional regression indicators:

1. **Before any code changes**, run a curated subset on current `main` and
   record baselines (final profiles, flame speed / consumption speed, heat
   release rate, run success):
   `example_single`, `example_diffusion`, `example_twin`,
   `example_cylindrical_outward`, `example_cylindrical_inward`,
   `example_laminarFlameSpeed`.
2. After Phase 1: rerun with `firstOrderUpwind` — results must match baselines
   to within run-to-run reproducibility (legacy path is bit-identical;
   residual differences only from thread scheduling).
3. Rerun with `secondOrderLimited` — expect small, physically consistent
   shifts; large or unphysical deviations are a red flag.

Baseline artifacts and comparison scripts live in `test/convergence/` (shared
with §6.4) with a README documenting how to regenerate them.

### 6.4 Grid-convergence studies (Phase 2)

Scripted and rerunnable in `test/convergence/`:

- **Case A:** unbounded strained premixed H₂ flame (FixedValue/ZeroGradient
  BCs).
- **Case B:** twin flame (ControlVolume left BC) and/or cylindrical flame
  (α=1, `rphalf` weighting) — exercises the remaining BC and geometry paths.

Protocol: sweep grid tolerances (`vtol`, `dvtol`, `gridMax`) to produce a
resolution ladder; compare consumption speed and peak temperature vs. N for
`firstOrderUpwind` vs. `secondOrderLimited`, against a fine-grid reference.
Include a trial with relaxed `dampConst` to quantify the achievable grid
coarsening. Record CVODE step counts for both schemes.

## 7. Phasing

- **Phase 0 — Baselines:** run and store example baselines (§6.3) from
  unmodified `main`. No code changes.
- **Phase 1 — Convection upgrade:** kernel + call-site integration + config
  plumbing + unit tests + regression parity + baseline comparison. Mergeable
  on its own.
- **Phase 2 — Convergence studies:** scripts + analysis (§6.4). Produces the
  **decision gate for Phase 3**: does diffusion error dominate total error at
  the target (coarsened) resolutions?
- **Phase 3 — Diffusion upgrade (conditional on Phase 2):**
  - *(a) minimal:* harmonic-mean face diffusivity within the existing 3-point
    tridiagonal structure; or
  - *(b) full:* 5-point 4th-order conservative flux stencils on the nonuniform
    grid + pentadiagonal extension of the BDF integrator (banded Thomas,
    still O(N)).
  Design details intentionally deferred until the gate is passed; only the
  decision criteria are fixed now.

## Appendix A: Future path — conservative flux-form framework (Approach C)

Out of scope for this effort, but on the radar as a concern with the current
formulation. Recorded here so later work can pick it up.

**Motivation.** The advective-form split operator (`∂y/∂t = −v ∂y/∂x + C`)
preserves constant fields trivially but is not discretely conservative:
species mass and enthalpy transport are not telescoping sums, so conservation
errors appear at the level of truncation error and interact with the
splitting constants. A face-flux (finite-volume) formulation would make
convection, diffusion, and the cross-term fluxes (jFick, jSoret, jCorr —
already face-based quantities in `updateCrossTerms`) share one machinery and
be conservative by construction.

**Sketch.**

1. **Staggered mass flux:** the continuity march already produces `rV` by
   integrating between nodes; reinterpret `rV` as a *face* quantity
   (`rV_{j+1/2}`), making the discrete continuity statement exact per cell:
   `(rV_{j+1/2} − rV_{j−1/2}) = −dlj_j · r_j · (drhodt + ρβU)_j`.
2. **Conservative update:** evolve `ρφ` (φ ∈ {Y_k, h, ...}) as
   `∂(ρφ)/∂t = −(1/r) d(rVφ̃)/dx + ...` with φ̃ the limited upwind face
   reconstruction from §3.1 — the same reconstruction machinery is reused.
3. **Constant preservation:** with fluxes `rV_{j+1/2}·φ̃_{j+1/2}` and the
   discrete continuity identity above, a constant φ field yields
   `∂(ρφ)/∂t = φ·∂ρ/∂t` exactly — the compatibility condition that makes
   conservative form safe under nonzero velocity divergence. This identity
   is the key invariant to unit-test.
4. **Complications to resolve when taken up:**
   - The `U` (momentum) equation has source structure (`−U²`, strain terms)
     and is naturally non-conservative; it can stay in advective form.
   - Splitting constants (`splitConst*`) are currently additive in `∂y/∂t`;
     in conservative variables they become additive in `∂(ρφ)/∂t`, changing
     how `SplitSolver` computes deltas and rebalancing terms.
   - Density is not a state variable of the convection sub-problem (it is
     reconstructed from `T`, `Wmx`); a conservative update of `ρφ` needs a
     consistent in-step model for `ρ(t)` — likely the same `drhodt` frozen
     field used by the continuity march.
   - Boundary closures (ControlVolume, WallFlux) are already flux balances
     and map naturally; ZeroGradient/outflow need one-sided flux closures.
5. **Suggested migration order:** species equations first (linear in φ, most
   to gain from conservation), then energy, leaving `U` advective. Reuse
   `ConvectionDifferencer`'s reconstruction internals by exposing the face
   values `ỹ_{j±1/2}` (a natural refactor seam already anticipated by the
   kernel design).

## Appendix B: Decisions log

| Decision | Choice | Rationale |
|---|---|---|
| Primary goal | Accuracy per grid point | Fewer points for same quality; enables coarser adaptive grids |
| Scope | Convection now; diffusion by experiment | Convection is 1st-order (dominant error); diffusion already 2nd-order |
| Scheme family | Limited 2nd-order upwind (MUSCL-type) | Bounded, robust, cheap stencil, standard for flame codes |
| Limiter | van Albada (default) | Smooth RHS for CVODE/Adams step-size control |
| Config | `convectionScheme` option; new default `secondOrderLimited`; legacy `firstOrderUpwind` retained bit-identical | A/B studies + rollback safety |
| Option naming | `secondOrderLimited` / `firstOrderUpwind` | States order up front; "Limited" flags local 1st-order clipping |
| Continuity march | Trapezoidal, tied to scheme option | 2nd-order; keeps legacy path bit-identical |
| Architecture | Shared kernel (Approach B) | DRY across UTW/species/quasi-2D call sites; unit-testable |
| Approach C | Deferred; documented in Appendix A | More than the current goal requires; concern noted for future |
| Validation | Grid convergence on strained + twin/cylindrical flames; example baselining | Covers all BC paths; compensates for incomplete unit suite |
| Quasi-2D | Bug fix accepted; behavior change acceptable | Orphan configuration, no tests, not a merge blocker |

## Phase 2 findings (Task 2.2 convergence study)

*Appended by the Task 2.2 analysis stage. Sections above are unchanged. All
numbers below are reproducible from the committed JSONs in
`test/convergence/results/` (main ladder) and
`test/convergence/results/diagnostics/` (anomaly probes); the matrix, the
failures, and the plots are documented in
`.superpowers/sdd/task-2.2-runs-report.md`. Runs were produced at commit
`f38cfc31`, `nThreads=1` (deterministic).*

**Method / reference convention.** Error is measured against the finest
*completed* `secondOrderLimited` (sol) run per case, matching
`plot_convergence.py` (strained: rung 4, N=223; twin: rung 5, N=525;
cylindrical: rung 3, N=187). This is a **self-reference**: sol's own error is
computed against itself, which understates sol's absolute error near the
reference N and inflates fitted sol orders. For strained and cylindrical the
sol ladder is also **truncated** (4 completed rungs each; §P2.4), so the sol
order estimates for those two cases rest on 4 and 3 points respectively and
should be read as "at least 2nd-order, converging fast," not as precise
orders. Two facts partly offset the self-reference worry for the scalar
*values*: (1) `firstOrderUpwind` (fou) approaches sol's converged value from
below in every case, and (2) the two schemes agree on the continuum limit to
≤0.05% on consumption_speed — so the finest-sol value is a sound truth proxy
even where its own error bar is unquantified.

### P2.1 Observed convergence orders

Log-log least-squares fit of relative error vs. N (order = −slope). "n pts"
notes how many ladder points entered the fit (reference run excluded).

| case | metric | sol order | fou order |
|---|---|---|---|
| strained | consumption_speed | 2.58 (4 pts) | 1.13 (6 pts) |
| strained | peak_T | *noise floor* (4 pts) | 0.94 (6 pts) |
| twin | consumption_speed | 3.07 (5 pts) | 0.95 (6 pts) |
| twin | peak_T | *noise floor* (5 pts) | 0.54 (6 pts) |
| cylindrical | consumption_speed | 3.18 (3 pts) | 1.55 (6 pts) |
| cylindrical | peak_T | *noise floor* (3 pts) | 1.06 (6 pts) |

- **fou behaves as ~1st-order** on both metrics in all cases (0.5–1.6),
  consistent with its O(h) numerical diffusion.
- **sol consumption_speed fits are ≥2.5**, i.e. *at least* 2nd-order. The
  absolute value is inflated by the self-reference and truncated ladders and
  should not be quoted as a precise order — the honest statement is "sol is
  at least 2nd-order and converges much faster than fou."
- **sol peak_T is at the noise floor** (relative error 1e-5…1e-4,
  non-monotonic) at *every* rung including the coarsest — it is grid-converged
  before the ladder even starts, so no order is extractable (the tabulated
  −1.37/1.65/1.93 are fits to noise and are meaningless).

### P2.2 Error magnitude and grid-size reduction at matched error

The order fits understate the real result; the per-point *error magnitude* is
the headline. Relative error at the coarsest sol grid vs. the best (finest)
fou grid:

| case | metric | sol @ coarsest | fou @ finest | grid reduction at matched error |
|---|---|---|---|---|
| strained | consumption_speed | 0.62% @ N=55 | 1.07% @ N=360 | **>6.5×** (fou's finest never reaches sol's coarsest accuracy) |
| strained | peak_T | 0.002% @ N=55 | 0.30% @ N=360 | **≫7×** (fou is ~100× less accurate even at finest) |
| twin | consumption_speed | 0.30% @ N=76 | 0.27% @ N=531 | **~7×** (fou needs ~N=500 to match sol's N=76) |
| twin | peak_T | 0.01% @ N=76 | 0.08% @ N=531 | **≫7×** (fou never matches within the ladder) |
| cylindrical | consumption_speed | 0.69% @ N=64 | 0.03% @ N=449 | **~1×** (fou's cs converges quickly here; no sol advantage) |
| cylindrical | peak_T | 0.05% @ N=64 | 0.12% @ N=449 | **>7×** (fou never matches within the ladder) |

**Reading:** the "accuracy per grid point" goal (§2) is strongly met for
**peak_T in all three cases** (sol is grid-converged at the coarsest grid;
fou never catches up within the tested ladder) and for **consumption_speed in
strained and twin** (~7× fewer points for matched error). The one honest
exception is **consumption_speed in the cylindrical case**, where fou's
consumption speed also converges quickly and sol shows no meaningful per-point
advantage on that single metric (its peak_T advantage there is still large).
So the win is real and consistent but not uniform across every (case, metric)
cell.

### P2.3 CVODE step-count / limiter-smoothness comparison (spec §5)

Total convection CVODE steps, sol/fou ratio along the ladder (matched N):

| case | rung0 | rung1 | rung2 | rung3 | rung4 | rung5 |
|---|---|---|---|---|---|---|
| strained | 1.11 | 1.17 | 1.23 | 1.66 | 1.88 | *sol failed* |
| twin | 0.88 | 1.06 | 1.14 | 1.24 | 2.36 | 2.76 |
| cylindrical | 1.06 | 1.12 | 1.28 | 1.29 | *sol failed* | *sol failed* |

The limited scheme's convection cost **grows super-linearly with refinement**
in every case — the sol/fou ratio rises monotonically from ~1 at the coarsest
grid. In `twin` (the only case where sol completes the whole ladder) it
plateaus at ~2.8× and the run is fine: an **elevated but survivable**
stiffness cost, directly attributable to the limiter's reduced RHS smoothness
(spec §5). In the two high-strain cases the cost does not plateau — it
diverges into the CVODE failure of §P2.4. Per-*global*-step convection counts
from the logs make the divergence concrete: twin sol at N=525 sits at a stable
~1,600 steps/global-step (vs. ~400 for fou); the failing strained sol rung 5
climbs 350→3,900 as its grid grows 350→505, and the failing cylindrical sol
rung 4 climbs into the 10,000–16,000 range, immediately before CVODE quits
(fou at comparable N stays at 70–720). **This is the spec §5 "step-count
blowup" that was pre-registered as a pre-merge red flag; it has triggered.**

### P2.4 Anomaly: deterministic CVODE failure at the finest rungs

**Facts (data shows).** `secondOrderLimited` fails deterministically with
`CVODE Integrator had too many errors` at strained rung 5 and cylindrical
rungs 4–5 (each failing every attempt across two invocations, `nThreads=1`).
`firstOrderUpwind` completes those same rungs at identical tolerances; `twin`
is unaffected for both schemes. In **every failing run** the adaptive grid
grows monotonically and never settles (strained rung 5: 167→505 across 23
regrids; cylindrical rung 4: 108→408) while the per-step convection CVODE
count climbs in lockstep, both diverging together in the final regrids before
the integrator gives up. In **every passing run** (both schemes, including
`twin` sol at N=525) the grid overshoots then plateaus at a stable point
count. Cells are nowhere near the `gridMin` floor (smallest ~2.3e-5 vs. floor
5e-7); refinement is driven by `vtol`/`dvtol` curvature demand around the
sharpening front.

**Diagnostic probes** (analysis-only, no `src/` edits; JSONs in
`results/diagnostics/`):

- **Grid frozen** (rung-5 tolerances, `regridStepInterval=1e9`,
  `regridTimeInterval=1e30`): **completes**, bounded convection cost (~42
  steps/global-step, N=100). The failure *requires the grid to keep
  refining* — a non-refining grid at the same tolerances integrates fine.
- **Intermediate "rung 4.5"** (`vtol=0.040, dvtol=0.065, gridMax=5.0e-5`):
  **completes** at N=269 with converged scalars (consumption_speed 0.35915,
  peak_T 1558.92) and a settled grid — a stable fixed point exists just below
  rung 5.
- (`Debug(regridding=False)` was found to be only a verbose-print toggle, not
  an adaptation switch, so that probe merely reproduced the failure
  identically — a useful determinism confirmation.)

**Hypotheses weighed.**

- *H1 — intrinsic limiter stiffness* (limiter non-smoothness stiffens the
  convection RHS): **supported but survivable.** twin's stable ~2.8× step-count
  elevation and sol>fou counts everywhere show it is real; it is the
  *mechanism* by which a run eventually exceeds CVODE's error budget, not by
  itself the trigger.
- *H2 — grid-adaptation feedback* (sol's sharper profiles keep triggering
  regrids): **best-supported as the proximate trigger.** The grid-frozen probe
  removes the failure; the monotonic non-settling grid growth is present in
  exactly the failing runs and absent in every passing run.
- *H3 — tolerance ladder simply too tight:* **refuted.** fou completes all
  rungs at identical tolerances with a settled grid, so the tolerances are not
  intrinsically unreasonable.

**Best-supported mechanism.** At the finest tolerances in high-strain flames
the low-numerical-diffusion limited scheme has **no stable grid fixed point**:
each refinement sharpens the resolved front (by design — that is the accuracy
gain), which both re-triggers refinement (interior insertion near the front;
boundary removal blocked by the `boundaryTolRm`=1e-5 flatness test) *and*
stiffens the convection sub-integration (H1), a positive feedback that
diverges until the per-step CVODE count exceeds its error budget and the
integrator quits. The threshold sits between "rung 4.5" (N=269, stable) and
rung 5 for strained. This is an **edge-of-stability phenomenon at the finest
one or two rungs, not a blanket breakdown** — and it is outside the *default*
operating envelope (default `vtol`=0.12 is coarser than rung 2; the failures
are at `vtol`=0.033–0.05, 2.4–3.6× tighter).

**What would resolve it (analysis-level, not a decision).** A durable fix is
code-level: (a) smooth the limiter near extrema and/or freeze σ within a
convection sub-step so the RHS CVODE sees is smoother (spec §5's own
hypothesis); or (b) a grid-adaptation guard that denies the runaway a foothold
(cap growth per regrid / add removal hysteresis). A definitive H1-vs-H2
isolation would need a run on a *fixed fine* grid (frozen at N≈500), which the
current harness cannot seed without a `src/` change.

### P2.5 dampConst trial (strained sol, dampConst 7 → 15)

Grid-point count is **identical point-for-point** (55/75/115/145/223 at rungs
0–4) between the default `dampConst=7` and the relaxed `dampConst=15`; rung 5
fails under both. **No additional coarsening headroom** from relaxing
`dampConst` for this case — the point count is governed by `vtol`/`dvtol`, not
the removal-damping term, at these strain/tolerance settings. The relaxed
value also did not rescue rung 5, consistent with §P2.4 (the failure is a
convection-integrator robustness issue, not a grid-oscillation/damping one).

### P2.6 Phase 3 gate assessment

*Gate criterion (from the plan): "does diffusion error dominate at target
resolutions? — if `secondOrderLimited` error-vs-N flattens toward
2nd-order-diffusion-limited behavior while convection refinement no longer
helps, the gate opens."*

- **The flattening signature is present.** sol's peak_T is grid-converged to
  ≤0.05% at the coarsest grid in all three cases, and sol's consumption_speed
  is within ~0.3–0.7% at the coarsest grid and flat thereafter. Refining the
  convection grid beyond N≈100–150 changes the sol scalars negligibly —
  convection refinement no longer helps.
- **Inference, with a caveat.** Because convection error is exhausted at
  coarsened target resolutions, the residual error there is *not*
  convection-limited; the 2nd-order diffusion operator is the leading
  remaining term, so diffusion error plausibly dominates. This is **inferred
  from convection-error exhaustion plus cross-scheme agreement, not measured**
  — the study did not A/B the diffusion scheme, so "diffusion dominates" is a
  well-motivated inference, not a direct measurement. The self-referenced
  metric also means the *absolute* error floor at target resolution is not
  quantified.
- **Independent robustness caveat.** The §5 step-count red flag has triggered
  (§P2.3–P2.4). This concerns the *convection scheme's* robustness under
  refinement, which is orthogonal to the diffusion-dominance gate criterion —
  it neither opens nor closes the gate, but it means Phase 2 is **not a clean
  bill of health** for the limited scheme at tight tolerances.

### P2.7 Recommendation (decision is the project owner's)

1. **Phase 3 gate: recommend OPEN on the diffusion-dominance criterion.**
   Convection error is exhausted at coarsened target resolutions and further
   convection refinement no longer helps, so a diffusion-accuracy upgrade is
   the next available lever. Proceed to the Phase 3 brainstorming/design pass
   (per the plan, a fresh scoped spec — do not improvise Phase 3 from the
   one-paragraph sketch).
2. **Condition the gate on resolving the §5 tight-tolerance CVODE robustness
   bug (§P2.4).** A convection scheme that *diverges* under grid refinement
   (rather than converging) is a latent robustness hazard even though it sits
   outside the default coarse-grid operating envelope, and it is exactly the
   red flag the spec pre-registered. Recommendation: fix it (option (a) or (b)
   above, chosen via a short design pass) **before Phase 3 diffusion work
   lands**, and — if the limited scheme ships as the default in the Phase 1
   merge — record it as a known limitation / release note in the meantime.
   This is a convection-scheme item, not a reason to keep the Phase 3 gate
   closed.
3. **dampConst: no action.** It yields no additional coarsening headroom for
   the strained case (§P2.5).

**Decision: PENDING user review.**
