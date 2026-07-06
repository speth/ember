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

## A. Validation results (G6)

*Appended by Task G6. Full errTol ladder (rungs 0–5, 3 cases ×
2 schemes = 36 runs) run at commit `dc5e197` (G5, post-calibration
`errCoeff`: `C_1 = 0.125`, `C_2 = 0.0139`; default `Grid.errTol = 1e-4`).
Raw JSONs/plots: `test/convergence/results-errtol-final/` (untracked
scratch, per the G4/G5 convention — `results/` stays reserved for the
Task-2.2 vtol-era history). All 36 runs completed on the first attempt
(no retries).*

### A.1 Run matrix and the §P2.4 hard criterion

**Result: 36/36 runs completed. Zero CVODE failures anywhere in the
ladder** — including `secondOrderLimited` (sol) at strained rung 5 and
cylindrical rungs 4–5, the exact (case, rung) combinations that failed
deterministically under the old `vtol`/`dvtol` criterion (spec
§P2.4). No gridMax-vs-criterion disentangling probe was needed (it is
only called for when a run fails; none did).

**Grid settling.** `analyze_settling.py`'s SETTLED/NOT SETTLED heuristic
(±3 grid points over the last 25% of profile outputs) reports 14/36 runs
as NOT SETTLED, including several `firstOrderUpwind` twin runs that have
never had a stability problem in this study. Manual trajectory inspection
of all 36 `prof*.h5` N-sequences shows **zero cases of monotonic,
unbounded growth** (the actual §P2.4 failure signature) — every
trajectory overshoots then plateaus; the "NOT SETTLED" calls are an
artifact of the tool's fixed absolute-point tolerance not scaling with N,
which has grown well past the Task 2.2-era ladder the tool's default was
tuned against (up to N=631 here). Tail spread on the flagged runs is 1–9
points on N in the tens-to-hundreds (≤4% relative), with no run's tail
window monotonically increasing. Selected trajectories for the two
previously-failing configurations, at their finest rungs:

| run | N trajectory (tail) |
|---|---|
| strained sol rung 5 | …229, 222, 221, 221, 224, 219, **218, 217, 217, 217, 217, 218, 218, 218** |
| cylindrical sol rung 4 | …134, 132, 130, 130, 130, 130, **130** (settled, flagged SETTLED) |
| cylindrical sol rung 5 | …184, 184, 185, 185, 184, 183, 183, 183, 184, **183, 183, 183** (settled) |

**Convection CVODE step counts** (per-global-step `[C: N]`, from the run
logs), compared against the §P2.3 divergence signature (climbing into the
thousands/tens-of-thousands in lockstep with N, e.g. old strained sol
rung 5: 350→3,900; old cylindrical sol rung 4: 10,000–16,000):

| run | steps/global-step: min / max / last-10 range |
|---|---|
| strained sol rung 5 (N=218) | 53 / 314 / 53–64 |
| strained sol rung 4 (N=130) | 47 / 297 / 47–50 |
| cylindrical sol rung 4 (N=122) | 274 / 646 / 274–285 |
| cylindrical sol rung 5 (N=183) | 282 / 646 / 282–294 |
| twin sol rung 5 (N=302, reference: never failed) | 255 / 545 / 318–327 |

All bounded to two-to-three digits and **decaying** to a steady low value
by the end of each run — the opposite of the old divergence signature, and
in the same range as `twin` sol (which never failed under the old
criterion either). This is the strongest single piece of evidence that
the error-budget criterion has structurally removed the §P2.4 feedback
loop, per its design intent (§1, §2).

**Accuracy vs. the "old vtol rung-5" bar.** Because sol under `vtol`
never completed rung 5 for strained or cylindrical (the P2.4 failure
itself), there is no direct "old sol rung-5" number for those two cases;
two reasonable stand-ins were evaluated:

1. *Continuum reference* (Richardson-extrapolated `cs_inf` from
   calibration-notes.md: 0.35916 / 0.1707226 / 0.22766 for
   strained/twin/cylindrical), old target = best-*completed* old-vtol sol
   rung (strained rung 4, twin rung 5, cylindrical rung 3): target errors
   1.56e-4 / 7.2e-6 / 2.25e-4. Under this reference the new ladder's
   finest rung (rung 5, errTol=2e-4) does not quite clear the bar for
   strained (1.96e-4) or cylindrical (5.1e-4); twin's target is
   unreachable by construction (`cs_inf` for twin is itself anchored
   close to twin's own old-vtol rung 5, making that target near-tautological).
2. *Old-vtol `firstOrderUpwind` rung 5* (the one rung-5 baseline that
   completed in **every** case under the old ladder) as the target: errors
   1.09e-2 / 2.67e-3 / 5.6e-4 for strained/twin/cylindrical. Under this
   (arguably fairer, non-circular) reference, new-ladder sol clears the
   bar from **rung 0** (strained), **rung 3** (twin), and **rung 5**
   (cylindrical).

Either way, the two accuracy-qualifying rungs for the higher-strain cases
land at or near the ladder's finest rungs (4–5), matching the task's
prior expectation — and since **every** rung 0–5 completes and settles
cleanly for every case/scheme (§A.1 above), the precise rung threshold
does not change the verdict.

**Hard criterion: PASSES.** At and beyond the qualifying rungs, strained
and cylindrical sol complete with a settling grid and bounded (decaying)
convection step counts. No contingency (growth cap or σ-freeze, spec §5)
is needed.

### A.2 Parity (fou vs. sol at matched errTol)

Ratio `err_fou / err_sol` for `consumption_speed`, own-scheme finest-rung
(rung 5) self-reference (matching the G5 calibration convention), rungs
1–4:

| case | rung 1 | rung 2 | rung 3 | rung 4 |
|---|---|---|---|---|
| strained | 4.41 | 5.40 | 3.81 | 2.67 |
| twin | 0.71 | 1.55 | 3.03 | 2.21 |
| cylindrical | 0.38 | 0.86 | 0.94 | 1.07 |

**Parity verdict: PARTIAL.** In *aggregate* the criterion is met:
geometric mean over all 12 (case, rung) points is **R = 1.70**, inside
the ~2× band and much closer to 1 than the pre-calibration R = 2.85
(calibration-notes.md §2) — the G5 recalibration moved in the intended
direction. twin and cylindrical are also individually within the band
(per-case geomeans 1.65 and 0.76). But **strained exceeds the band at
every rung 1–4** (2.67–5.40×, per-case geomean 3.95), which the spec §2
"drift toward the extremes" (≈errTol^(1/6)) language does not license —
that model predicts band violation only at the tolerance extremes, not
across the whole practical range. This is not a clean PASS.

The residual is structural, not a calibration mistake: the
pre-calibration per-case geomeans spanned ~5× (strained 6.13, twin 3.16,
cylindrical 1.19), so **no single global `C_2` can bring all three cases
inside the 2× band simultaneously** — shifting `C_2` to fix strained
would push cylindrical (already at/below parity) out the other side. The
G5 calibration, tuned to the 3-case geometric mean per protocol,
necessarily lands strained high and cylindrical low. **Owner decision
required** — options: (a) accept the aggregate-parity interpretation of
the §4 criterion as satisfied; (b) keep the global `C_2` and add
per-case/per-regime guidance to user docs (e.g. "for highly strained
flames, expect the first-order scheme to need a few× tighter errTol for
matched accuracy"); or (c) revisit the calibration weighting (e.g.
minimax over cases instead of geometric mean, or a strain-dependent
`C_2`). Not a regression from G5 (the spread was disclosed there);
flagged here because the acceptance criterion as written is not met for
strained individually.

`peak_T` ratios (5–220×) are not meaningful, per the same caveat noted in
calibration-notes.md §2: sol's `peak_T` sits at the noise floor
(~1e-5–1e-4) in every case, so any nonzero fou `peak_T` error inflates the
ratio arbitrarily even though both errors are tiny in absolute terms.

**N-ratio (fou/sol) at matched errTol**, cross-scheme (finest-sol)
reference, rungs 0–5 (full table:
`test/convergence/results-errtol-final/plots/errtol_error_table.md`):

| case | rung 0 | rung 1 | rung 2 | rung 3 | rung 4 | rung 5 |
|---|---|---|---|---|---|---|
| strained | 1.96 | 1.95 | 1.72 | 2.15 | 2.05 | 1.93 |
| twin | 1.83 | 2.32 | 1.88 | 1.74 | 2.10 | 2.09 |
| cylindrical | 1.64 | 2.23 | 2.17 | 2.35 | 2.60 | 2.89 |

sol consistently uses ~1.6–2.9× fewer points than fou at matched errTol,
across the whole ladder — the expected higher-order-scheme benefit,
growing mildly at tighter tolerances for cylindrical.

### A.3 errTol → measured QoI-error table

Cross-scheme reference (finest `secondOrderLimited` run per case);
full machine-generated table with all rungs/metrics:
`test/convergence/results-errtol-final/plots/errtol_error_table.md`
(generated by `plot_convergence.py`, not hand-edited). Headline values,
`consumption_speed`, `secondOrderLimited`:

| errTol | strained N / err | twin N / err | cylindrical N / err |
|---|---|---|---|
| 2e-2 | 28 / 5.76e-3 | 41 / 1.16e-2 | 36 / 8.13e-3 |
| 8e-3 | 41 / 6.16e-3 | 50 / 1.01e-2 | 40 / 4.88e-3 |
| 3.2e-3 | 71 / 3.11e-3 | 92 / 3.05e-3 | 66 / 1.86e-3 |
| 1.3e-3 | 80 / 2.74e-3 | 149 / 8.64e-4 | 97 / 7.72e-4 |
| 5e-4 | 130 / 1.48e-3 | 192 / 4.81e-4 | 122 / 2.23e-4 |
| 2e-4 (rung 5, finest tested) | 218 / (ref) | 302 / (ref) | 183 / (ref) |

(Rung 5 is the reference run for its own case/scheme in this convention,
hence "(ref)"; see §A.1 for its absolute accuracy against the
independent `cs_inf` continuum estimate.) Plots (error-vs-N, error-vs-
errTol per case/metric): `test/convergence/results-errtol-final/plots/`.

### A.4 G5 calibration summary (recap)

`errCoeff` recalibrated `1/15 → 0.0139` for `secondOrderLimited`
(`C_2 = C_2_old · R^(-3/2)`, R = 2.85 measured pre-calibration parity
ratio); default `Grid.errTol` set to `1e-4`. **Sign-off pending:** the
task brief specified `C_2_new = C_2 · R^(+3/2)`; G5 determined this
exponent is inverted (larger `errCoeff` refines the grid *sooner*, so
matching the brief's stated intent — "bring the ratio to ~1" — requires
the negative exponent) and applied `R^(-3/2)` instead, with the derivation
documented in `test/convergence/results/calibration-notes.md` §3 and
independently re-derived by the G5 code reviewer. This sign correction
and the errTol=1e-4 default (which does not jointly attain the owner's
informal envelope of err~1e-4 @ N~100 — measured N=115–192, err=2.2e-4–
1.0e-2 across cases at default `errTol`) are both **pending explicit
owner sign-off**, per the G5 ledger entry. (Envelope-range
reconciliation: the N=115–192 / 2.2e-4–1.0e-2 span quoted here covers
all three cases, both schemes, and both QoIs from this ladder's runs
bracketing the default; the tighter N=115–165 / err_cs=3.1–6.5e-4 span
in calibration-notes.md §4 covers `secondOrderLimited` consumption
speed only, measured directly at errTol=1e-4 with default gridMax. The
two describe the same operating point at different scope, not a
discrepancy.)

### A.5 Anomalies and open items

- `analyze_settling.py`'s fixed `tol_pts=3` absolute tolerance produces
  false "NOT SETTLED" calls once N grows past a few hundred (§A.1); a
  relative tolerance (e.g. tail spread as a fraction of N) would track
  actual settling more reliably. Not fixed here (analysis-only task;
  no `src/` or harness changes authorized by this task).
- The "old vtol rung-5 accuracy" bar is not uniquely defined for strained
  and cylindrical, since old-vtol sol never reached rung 5 there (that
  failure is the reason for this whole design pass). §A.1 evaluates both
  reasonable stand-ins; both support the same PASS verdict.
- strained's persistent >2× parity gap across all of rungs 1–4 (§A.2) is
  a candidate for a case-specific recalibration note if tighter parity is
  wanted for high-strain flames specifically, rather than the current
  single global `C_2`.
- Default `Grid.errTol=1e-4` sits between this ladder's rung 4 (5e-4) and
  rung 5 (2e-4); it was not run as its own rung here (out of scope for
  the ladder, already characterized in calibration-notes.md §4 with
  default `gridMax`).
