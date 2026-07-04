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
