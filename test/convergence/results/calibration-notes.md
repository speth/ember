# Task G5: Calibration of C_2 (secondOrderLimited errCoeff) and the default errTol

Working notes for the calibration pass. Raw JSONs for the runs referenced here
live in `test/convergence/results-errtol/` (pre-calibration, C_2 = 1/15,
rungs 1-4 of the errTol ladder, all 3 cases x both schemes) and
`test/convergence/results-errtol-recal/` (post-calibration runs, C_2 = 0.0139).
These directories are untracked scratch outputs; this file records the numbers
that matter. All runs completed on the first attempt (no retries) unless noted.

Code state: errTol-ladder harness from G4; error-budget criterion from G2
(`E = errCoeff * h^(p+1) * max|d^(p+1)v| <= errTol * range(v)`), with
firstOrderUpwind (fou): p=1, C_1 = 0.125 and secondOrderLimited (sol): p=2,
C_2 = 1/15 before this calibration.

## 1. Pre-calibration convergence tables (C_2 = 1/15)

Reference for the error columns = rung 4 of the same (case, scheme), per the
Phase-2 protocol. `err_cs` / `err_pT` are relative errors of
`consumption_speed` / `peak_T`.

### strained / firstOrderUpwind

| rung | errTol | N | consumption_speed | peak_T | err_cs | err_pT |
|------|--------|---|-------------------|--------|--------|--------|
| 1 | 8e-3   | 80  | 0.3469285 | 1543.80 | 2.33e-02 | 6.48e-03 |
| 2 | 3.2e-3 | 122 | 0.3506173 | 1548.81 | 1.29e-02 | 3.26e-03 |
| 3 | 1.3e-3 | 172 | 0.3528842 | 1550.86 | 6.51e-03 | 1.94e-03 |
| 4 | 5e-4   | 266 | 0.3551978 | 1553.87 | (ref) | (ref) |

### strained / secondOrderLimited

| rung | errTol | N | consumption_speed | peak_T | err_cs | err_pT |
|------|--------|---|-------------------|--------|--------|--------|
| 1 | 8e-3   | 58  | 0.3577523 | 1559.53 | 3.73e-03 | 2.65e-04 |
| 2 | 3.2e-3 | 87  | 0.3580337 | 1559.20 | 2.94e-03 | 5.20e-05 |
| 3 | 1.3e-3 | 106 | 0.3588126 | 1559.19 | 7.74e-04 | 4.39e-05 |
| 4 | 5e-4   | 156 | 0.3590904 | 1559.12 | (ref) | (ref) |

### twin / firstOrderUpwind

| rung | errTol | N | consumption_speed | peak_T | err_cs | err_pT |
|------|--------|---|-------------------|--------|--------|--------|
| 1 | 8e-3   | 116 | 0.1723225 | 1839.22 | 6.07e-03 | 9.72e-04 |
| 2 | 3.2e-3 | 173 | 0.1719088 | 1840.32 | 3.66e-03 | 3.75e-04 |
| 3 | 1.3e-3 | 259 | 0.1715493 | 1840.79 | 1.56e-03 | 1.16e-04 |
| 4 | 5e-4   | 403 | 0.1712826 | 1841.01 | (ref) | (ref) |

### twin / secondOrderLimited

| rung | errTol | N | consumption_speed | peak_T | err_cs | err_pT |
|------|--------|---|-------------------|--------|--------|--------|
| 1 | 8e-3   | 83  | 0.1702724 | 1842.79 | 2.47e-03 | 1.90e-06 |
| 2 | 3.2e-3 | 116 | 0.1704729 | 1842.65 | 1.30e-03 | 7.29e-05 |
| 3 | 1.3e-3 | 183 | 0.1706369 | 1842.78 | 3.38e-04 | 5.12e-06 |
| 4 | 5e-4   | 234 | 0.1706947 | 1842.79 | (ref) | (ref) |

### cylindrical / firstOrderUpwind

| rung | errTol | N | consumption_speed | peak_T | err_cs | err_pT |
|------|--------|---|-------------------|--------|--------|--------|
| 1 | 8e-3   | 89  | 0.2271144 | 1781.83 | 1.61e-03 | 4.01e-03 |
| 2 | 3.2e-3 | 143 | 0.2271694 | 1785.27 | 1.36e-03 | 2.09e-03 |
| 3 | 1.3e-3 | 228 | 0.2273687 | 1787.73 | 4.89e-04 | 7.13e-04 |
| 4 | 5e-4   | 317 | 0.2274799 | 1789.00 | (ref) | (ref) |

### cylindrical / secondOrderLimited

| rung | errTol | N | consumption_speed | peak_T | err_cs | err_pT |
|------|--------|---|-------------------|--------|--------|--------|
| 1 | 8e-3   | 61  | 0.2271762 | 1791.19 | 1.91e-03 | 3.64e-04 |
| 2 | 3.2e-3 | 93  | 0.2274087 | 1791.61 | 8.90e-04 | 1.30e-04 |
| 3 | 1.3e-3 | 126 | 0.2275270 | 1791.80 | 3.71e-04 | 2.04e-05 |
| 4 | 5e-4   | 175 | 0.2276114 | 1791.84 | (ref) | (ref) |

## 2. Parity ratios and R

`ratio(rung) = err_fou / err_sol` at matched errTol, QoI = consumption speed
(rung-4 same-scheme references; rung 4 is excluded since its self-error is 0):

| case | rung 1 | rung 2 | rung 3 |
|------|--------|--------|--------|
| strained    | 6.25 | 4.38 | 8.42 |
| twin        | 2.45 | 2.81 | 4.60 |
| cylindrical | 0.84 | 1.53 | 1.32 |

**R = geometric mean over 9 (case, rung) pairs = 2.85** — outside the
[0.5, 2] keep-band, so C_2 was adjusted (section 3).

Sanity checks:

- peak_T: R_pT = 30.8, same direction (fou worse everywhere) but the
  magnitude is not meaningful — sol peak_T errors run down to ~2e-6, i.e.
  at reference-noise level, which inflates the ratio arbitrarily.
- Reference contamination: the rung-4 same-scheme reference biases the
  slow-converging fou errors low much more than the sol errors. Recomputing
  the ratios against scheme-independent continuum limits (Richardson
  extrapolation of the prior vtol-study sol ladders: cs_inf = 0.35916 /
  0.1707226 / 0.22766 for strained / twin / cylindrical) gives R = 4.61.
  The calibration below follows the specified protocol value R = 2.85, so it
  is conservative: after recalibration sol remains somewhat *more* accurate
  than fou at matched errTol on average (and cylindrical, which was already
  at parity, is not pushed far below it).

## 3. C_2 decision

Applied: **C_2 = 1/15 -> 0.0139** in `OneDimGrid::setOptions`
(src/grid.cpp), i.e. C_2_new = (1/15) * R^(-3/2) = (1/15)/4.81.

Note on the decision rule as written in the task brief: the brief stated
`err_sol ~ C_2^(2/3)` and `C_2_new = C_2 * R^(3/2)`. The proportionality has
an inverted sign: in the adaptation criterion a *larger* errCoeff produces a
*finer* grid (E exceeds the budget sooner), so `err_sol ~ (errTol/C_2)^(2/3)`.
Multiplying the sol error by R (the brief's stated intent, "bringing the
ratio to ~1") therefore requires scaling C_2 by R^(-3/2), not R^(+3/2).
Applying the +3/2 exponent literally would have refined sol grids ~1.7x
further and pushed the measured ratio from 2.85 up to ~8, away from parity.
The confirmation runs below verify the applied direction moves the ratio as
intended.

GridAdaptation gtests after the change: all 7 pass
(`./bin/unittest --gtest_filter='GridAdaptation.*'`).

### Confirmation runs (strained, rungs 2-3, both schemes, C_2 = 0.0139)

fou control: bit-identical to pre-calibration (C_1 untouched) — confirmed
for rungs 2 and 3 (same N, same consumption speed to all digits).

sol, errors vs the continuum limit cs_inf = 0.35916 (the post-cal rung-4 run
is too coarse, N = 130, to serve as a reference for its own ladder):

| rung | N pre -> post | err_cs pre -> post | factor (predicted x2.85) |
|------|---------------|--------------------|--------------------------|
| 2 | 87 -> 71   | 3.14e-3 -> 2.92e-3 | x0.93 |
| 3 | 106 -> 80  | 9.67e-4 -> 2.54e-3 | x2.63 |
| 4 | 156 -> 130 | 1.94e-4 -> 1.28e-3 | x6.62 |

Grids coarsen at every rung as expected. The error moves are noisy per rung
(the strained pre-cal curve is lumpy: rung 2 was barely better than rung 1,
rung 4 anomalously good), but the geometric-mean move over rungs 2-4 is
x2.5, matching the predicted x2.85 within the scatter, and the matched-errTol
ratio at rung 3 drops from 18.1 to 6.9 (continuum refs). Direction and
magnitude confirmed.

## 4. Default errTol

Decision runs: sol scheme, post-calibration C_2 = 0.0139, **default
gridMax = 2e-4** (unlike the ladder, which co-scales gridMax — this is what
a default-config user actually gets). Errors vs the continuum limits given
above (their own uncertainty is ~1.7e-4 relative for strained, ~1.2e-5 for
twin, ~1.8e-4 for cylindrical).

| errTol | strained N / err_cs | twin N / err_cs | cylindrical N / err_cs | geo-mean err |
|--------|--------------------|-----------------|------------------------|--------------|
| 2e-4 | 98 / 1.12e-3  | 141 / 7.30e-4 | 106 / 6.14e-4 | 7.9e-4 |
| **1e-4** | **115 / 6.5e-4** | **165 / 4.2e-4** | **154 / 3.1e-4** | **4.4e-4** |
| 5e-5 | 146 / 4.11e-4 | 205 / 3.29e-4 | 172 / 1.51e-4 | 2.7e-4 |
| 2e-5 | 184 / 1.51e-4 | (not run) | (not run) | — |

**Chosen default: errTol = 1e-4** (`Grid.errTol` in python/ember/input.py;
docstring examples: high accuracy 2e-5, minimal accuracy 5e-4, i.e. 5x
tighter / 5x looser).

Rationale and honest accounting against the owner envelope
(err_cs ~ 1e-4 at N ~ 100):

- The envelope is not jointly attainable on these cases with this scheme:
  reaching err_cs ~ 1e-4 measurably requires errTol ~ 2e-5 and N ~ 180-250,
  while N ~ 100 (errTol = 2e-4) delivers err_cs ~ 8e-4.
- errTol = 1e-4 is the balanced one-significant-figure point: N stays within
  1.2-1.7x of 100 while errors land at 3-6.5e-4 (part of which is reference
  uncertainty for strained/cylindrical). Tightening to 5e-5 buys only ~1.6x
  in error for ~25% more points, with visible flattening on twin
  (4.2e-4 -> 3.3e-4), consistent with a far-field error floor from the
  default gridMax.
- Robustness: with the ladder's co-scaled gridMax = 4e-5, strained fails
  deterministically (CVODE "too many errors", 3/3 attempts) at
  errTol = 2.5e-5 and 5e-5 — the known stability boundary in that regime.
  With the *default* gridMax = 2e-4 the same case runs cleanly down to at
  least errTol = 2e-5 (N = 184, 2.5 s). The default 1e-4 sits well inside
  the stable region; the docstring's high-accuracy example 2e-5 was verified
  to run on the strained case with default gridMax.

## 5. Files changed

- `src/grid.cpp`: sol-branch `errCoeff` 1/15 -> 0.0139 (+ comment).
- `python/ember/input.py`: `Grid.errTol` default 2e-3 -> 1e-4; docstring
  updated (measured accuracy statement; example values 2e-5 / 5e-4).
