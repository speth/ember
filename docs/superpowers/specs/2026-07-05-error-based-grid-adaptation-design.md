# Error-based grid adaptation (equidistribution design pass)

**Date:** 2026-07-05
**Branch:** `convection-scheme`
**Status:** Design approved by project owner; supersedes the vtol/dvtol
adaptation criteria as the robustness/Phase-3-precondition design pass
called for in the convection spec addendum §P2.7 and the SDD ledger
scoping note.

**Amended 2026-07-05 (pre-plan):** (1) error-estimate exponent corrected
from `h^p` to `h^(p+1)` — the approved `h^p` form has units of v/length
and cannot be compared to `errTol·range(v)`; the corrected form is the
standard degree-p interpolation-error bound and strengthens the
fixed-point argument. (2) A minimal fix for a pre-existing bug is in
scope: `adapt()`'s damping criterion reads `dampVal` after
`updateValues()` resizes it (Eigen resize does not preserve contents),
i.e. uninitialized values once a point has been inserted. (3) The
regression criterion is restated in accuracy-envelope terms (the ~1e-6
reproducibility floor only applies to identical configurations, which a
default-tolerance change is not).

## 1. Problem and goals

Two problems, one root cause:

1. **§P2.4 CVODE runaway.** At tight `vtol`/`dvtol`, `secondOrderLimited`
   fails deterministically (strained finest rung, cylindrical finest two):
   the grid refines monotonically without settling while convection CVODE
   step counts diverge. Diagnosed mechanism: the relative-change tests are
   *scale-free* — `|v[j+1]-v[j]| > vtol·range(v)` compares adjacent-point
   jumps, and `range(dv)` itself grows as the resolved front sharpens — so
   refinement can re-trigger indefinitely. `firstOrderUpwind` masks this
   with numerical diffusion.
2. **Tolerances don't mean accuracy.** `vtol`/`dvtol` are resolution
   heuristics: the same values produce very different solution accuracy
   for the two convection schemes, so users comparing schemes (or carrying
   old tolerance habits forward) draw wrong conclusions and over-refine.

**Goals:**

- A single user tolerance with *accuracy* semantics: same value → similar
  QoI accuracy under either convection scheme, with the higher-order
  scheme using fewer points. This makes the new scheme's benefit directly
  visible and prevents old-tolerance carry-forward.
- Structurally remove the §P2.4 runaway: an error-budget criterion has a
  grid fixed point by construction (per-cell error shrinks as h^(p+1)
  while the derivative estimate converges to the true profile's finite
  value).

**Owner decisions incorporated:**

- Tolerance semantics: **range-normalized relative** local error (one
  dimensionless knob; no per-species atol; near-zero minor species
  handled by the existing `absvtol` floor).
- Migration: **replace with deprecation** — new criterion is the only
  documented path; `vtol`/`dvtol` parse with a warning and are ignored.
- Robustness: fixing the §P2.4 failures is a **hard acceptance
  criterion**; limiter σ-freeze/smoothing is in scope as a contingency if
  grid changes alone don't cure them.
- Full equidistribution remesh and adjoint/QoI refinement rejected /
  out of scope (QoI work earmarked for Cantera).

## 2. Error estimator and adaptation criterion

For each adapted variable `v_k` and grid interval `j`:

```
E[k][j] = C_p · h_j^(p+1) · W[k][j]
```

- `p` = convection scheme spatial order: 1 (`firstOrderUpwind`),
  2 (`secondOrderLimited`).
- `W[k][j]` estimates `|d^(p+1) v_k / dx^(p+1)|` near interval `j`,
  built by repeated application of the existing nonuniform
  first-derivative stencils (`cfp/cf/cfm`, already computed in
  `adaptGrid`): `dv` as today, divided differences of `dv` for the
  second derivative, once more for the third when `p = 2`.
- **Max-filter over ±1 neighboring intervals** applied to `W` to control
  divided-difference noise near sharp fronts; biases conservatively
  toward refinement and prevents monitor flicker. At domain ends, the
  nearest interior value is used where a stencil application loses a
  point.

**Insertion** (replaces the two tests at `grid.cpp:123` and
`grid.cpp:133`):

```
insert if  E[k][j] > errTol · range(v_k)
```

The `range(v_k) < absvtol` minor-species skip is unchanged. The dvtol
gradient test is *subsumed*, not ported — the (p+1)-th derivative
measures what it stood in for.

**Removal** (same hysteresis structure as today): remove point `j` only
if the merged interval's error `E_merged[k][j]` (evaluated with
`h = hh[j]+hh[j-1]`) is below `rmTol · errTol · range(v_k)` for all `k`.
Merging ≈doubles `h`, so `E_merged ≈ 2^(p+1) ×` the sub-interval
errors — hysteresis is inherently stronger than the old criterion;
`rmTol` keeps its current role and default.

**Calibration.** `C_p` start at the classical interpolation-error
bounds — `C_1 = 1/8` (piecewise-linear), `C_2 = 1/15` (≈ piecewise-
quadratic) — and receive one empirical calibration pass on the three
study cases (strained/twin/cylindrical), chosen so matched `errTol`
yields matched measured QoI error (consumption speed, peak T) across
schemes at the practical operating point. Exact parity at *every*
tolerance is not achievable with schemes of different order: with this
monitor, QoI error scales ≈ `errTol^(p/(p+1))`, so the cross-scheme
error ratio drifts slowly (≈ `errTol^(1/6)`) away from the calibration
point. The measured `errTol` → QoI-error table (§4) documents the
actual mapping; the ~2× parity band in §4 covers the drift over the
practical tolerance range.

**Stated assumption.** For `p = 2` the monitor uses `h³·|v‴|` (the
quadratic-representation error associated with the convection scheme)
as a proxy for total local error; no separate fourth-derivative term is
included for the 2nd-order diffusion operator's error, since
fourth-derivative estimation from grid data is too noisy to be worth it
and `C_2` absorbs the difference. `p` stays scheme-based regardless of
any Phase-3 diffusion-order change.

**Untouched machinery:** damping (`dampVal`/`dampConst`), uniformity
tests, `gridMin`/`gridMax`, center/boundary special cases, boundary
add/remove logic, and the overall insert-then-remove pass structure.
One exception (pre-existing bug, minimal in-scope fix): `adapt()` must
maintain a working copy of `dampVal` that is kept consistent as points
are inserted/removed, because `updateValues()` resizes `dampVal`
without preserving its contents, leaving the damping criterion reading
uninitialized values mid-pass once the grid has changed.

## 3. Config surface and plumbing

- **New key `grid.errTol`** — the single range-normalized local-error
  tolerance. Default chosen during calibration so `secondOrderLimited`
  at defaults lands in the owner's accuracy envelope (scalars to ~1e-4,
  N≈100 on the study cases).
- **Deprecation:** `vtol`/`dvtol` still parse; setting either emits a
  warning (Python side, `input.py` validation, naming `errTol` as the
  replacement) and the values are ignored. Docstrings and examples
  updated to `errTol` only. `absvtol` unchanged.
- **Plumbing:** `readConfig` gains `errTol`; the grid error order is
  derived from the existing `convectionScheme` option and passed into
  `oneDimGrid` alongside the other adaptation parameters (same pathway
  as `vtol_in`, `grid.cpp:18`).
- **Input validation:** `errTol > 0`; warn above ~0.5 (criterion
  degenerates to "never refine").

## 4. Testing and acceptance

**Unit tests (gtest, alongside the Phase-1 convergence suite; TDD):**

1. Divided-difference estimator converges at the expected rate on
   analytic profiles (tanh) on a *nonuniform* grid; max-filter behaves
   at domain ends.
2. `adaptGrid` on a synthetic tanh profile meets the budget: post-adapt
   `max_j E/range ≤ errTol`.
3. No-flicker/idempotence: a second adapt pass changes nothing; removed
   points are not immediately re-inserted.
4. Order scaling: same profile and `errTol`, p=2 yields materially fewer
   points than p=1.

**Acceptance validation (Task 2.1 harness; strained/twin/cylindrical;
both schemes; `errTol` ladder):**

- **Hard criterion:** at `errTol` values reaching or exceeding
  old-rung-5 accuracy, all runs complete with a *settling* grid
  (overshoot-then-plateau, no monotonic growth) — the §P2.4 fix
  demonstrated.
- **Parity:** at matched `errTol`, fou and sol QoI errors agree within
  ~2×; sol/fou point-count ratio reported.
- **Docs artifact:** measured `errTol` → QoI-error table.
- **Regression:** default-settings example runs reproduce the phase-1
  baseline QoI scalars within the owner's accuracy envelope (consumption
  speed and peak T within 0.5%); N same or fewer. (The ~1e-6
  reproducibility floor applies only to identical configurations and is
  not the standard here, since the default tolerance necessarily
  changes.)

## 5. Contingencies (pre-scoped, not pre-built)

- **Per-regrid growth cap** — only if validation shows residual grid
  creep despite the fixed-point argument.
- **Limiter σ-freeze / smoothing** — only if a previously failing
  configuration still fails *with a settled grid* (pure H1 stiffness);
  joins the plan as a scoped item, not a new design pass.
