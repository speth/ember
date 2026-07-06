# Error-Based Grid Adaptation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace Ember's vtol/dvtol relative-change grid adaptation with a scheme-aware local-error budget (`grid.errTol`), fixing the §P2.4 tight-tolerance CVODE runaway and making the same tolerance mean the same accuracy under both convection schemes.

**Architecture:** A new per-interval error estimate `E = C_p · h^(p+1) · |d^(p+1)v|` (p = convection scheme order) replaces the two relative-change tests inside the existing `OneDimGrid::adapt()` insert/remove framework. Everything else in adaptation (damping, uniformity, min/max spacing, boundary logic, `rmTol` hysteresis structure) is preserved, except a minimal fix so `dampVal` stays consistent while the grid changes mid-pass. Python config gains `errTol`; `vtol`/`dvtol` are deprecated and ignored.

**Tech Stack:** C++ (Eigen, GoogleTest via `scons`), Cython bridge (`_ember.pxd`/`_ember.pyx`), Python config layer (`input.py`, pytest), convergence harness in `test/convergence/`.

**Spec:** `docs/superpowers/specs/2026-07-05-error-based-grid-adaptation-design.md` (as amended 2026-07-05 — note the `h^(p+1)` exponent).

## Global Constraints

- Branch: `convection-scheme`. **No push, no PR** (owner decision; branch is held locally).
- Both convection schemes (`firstOrderUpwind`, `secondOrderLimited`) must remain fully working; the error order follows the scheme.
- Do not touch the diffusion operator (that is Phase 3, separately gated).
- `absvtol`, `rmTol`, damping, uniformity, `gridMin`/`gridMax`, center/boundary logic keep their current roles; only the two vtol/dvtol tests are replaced (plus the scoped `dampVal` consistency fix).
- Initial estimator constants: `C_1 = 1/8` (p=1), `C_2 = 1/15` (p=2); Task G5 calibrates `C_2` and the default `errTol`. Provisional `errTol` default: `2e-3`.
- Commit messages: `convection: [G<N>] <summary>`.
- Build: `pixi run -- scons bin/unittest` (unit-test binary only) or `pixi run test` (full C++ + Python suite). If running `bin/unittest` directly fails with library-path errors, fall back to `pixi run -- scons test`.
- The hard acceptance criterion (Task G6): configurations at accuracy targets at/beyond the old failing rungs (strained rung 5, cylindrical rungs 4–5) must complete with a settling grid. If they don't, STOP and report — contingencies (growth cap, limiter σ-freeze) are an owner decision per spec §5.

## File Structure

- `src/grid.h` / `src/grid.cpp` — `errTol`/`errorOrder`/`errCoeff` members, `computeErrorWeights()`, new insert/remove tests, `dampVal` consistency fix. (Removes `vtol_in`, `dvtol_in`, `vtol`, `dvtol`.)
- `src/readConfig.h` — `errTol` field replaces `vtol`/`dvtol`.
- `python/ember/_ember.pxd`, `python/ember/_ember.pyx` — struct field + option transfer.
- `python/ember/input.py` — `Grid.errTol` option, vtol/dvtol deprecation, `Config._warnDeprecated()`.
- `python/ember/examples/example_laminarFlameSpeed.py` — only example that sets vtol/dvtol.
- `test/test_gridAdaptation.cpp` — new gtest file (picked up automatically by the `Glob('build/test/*.cpp')` in SConstruct).
- `test/python/test_config.py` — deprecation-warning tests.
- `test/convergence/run_convergence.py`, `test/convergence/plot_convergence.py`, new `test/convergence/analyze_settling.py` — errTol ladder, settling analysis.

---

### Task G1: Error-weight estimator (`computeErrorWeights`)

**Files:**
- Modify: `src/grid.h` (members near line 48, method declaration near line 185)
- Modify: `src/grid.cpp` (implementation after `setOptions`)
- Create: `test/test_gridAdaptation.cpp`

**Interfaces:**
- Produces: `void OneDimGrid::computeErrorWeights(const dvector& v, dvector& W) const` — fills `W` (size `jj+1`) with nodal estimates of `|d^(errorOrder+1) v / dx^(errorOrder+1)|`; requires `updateValues()` to have been called. New public members `double errTol; int errorOrder; double errCoeff;` (used by G2's criterion). Existing `vtol_in`/`dvtol_in`/`vtol`/`dvtol` members are NOT removed yet (G2 does that) so the build stays green.

- [ ] **Step 1: Write the failing test**

Create `test/test_gridAdaptation.cpp`:

```cpp
#include "../src/grid.h"
#include "../src/mathUtils.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>

namespace {

const double kXcenter = 1e-3;  // tanh front center [m]
const double kDelta = 1e-4;    // tanh front thickness [m]

double tanhProfile(double xx) {
    return std::tanh((xx - kXcenter) / kDelta);
}

//! Planar grid with FixedValue BCs and adaptation parameters chosen so that
//! only the criteria under test can fire (damping/gridMax neutralized).
OneDimGrid makeAdaptGrid(const dvec& x)
{
    OneDimGrid grid;
    grid.alpha = 0;
    grid.beta = 1;
    grid.leftBC = BoundaryCondition::FixedValue;
    grid.rightBC = BoundaryCondition::FixedValue;
    grid.setSize(x.size());
    grid.x = x;
    grid.absvtol = 1e-10;
    grid.rmTol = 0.6;
    grid.uniformityTol = 2.5;
    grid.gridMin = 1e-8;
    grid.gridMax = 1e-1;
    grid.dampConst = 7;
    grid.centerGridMin = 1e-8;
    grid.fixedLeftLoc = false;
    grid.errTol = 2e-3;
    grid.errorOrder = 2;
    grid.errCoeff = 1.0 / 15.0;
    grid.updateValues();
    grid.dampVal = dvec::Constant(grid.nPoints, 1e6);  // damping never binds
    return grid;
}

dvec uniformGrid(size_t n, double x0, double x1)
{
    dvec x(n);
    for (size_t j = 0; j < n; j++) {
        x[j] = x0 + (x1 - x0) * j / (n - 1);
    }
    return x;
}

dvector sampleProfile(const OneDimGrid& grid)
{
    dvector v(grid.nPoints);
    for (size_t j = 0; j < grid.nPoints; j++) {
        v[j] = tanhProfile(grid.x[j]);
    }
    return v;
}

double maxWeight(OneDimGrid& grid, const dvector& v)
{
    dvector W;
    grid.computeErrorWeights(v, W);
    return *std::max_element(W.begin(), W.end());
}

} // namespace

// max|v''| = (4/(3*sqrt(3)))/delta^2 for tanh; max|v'''| = 2/delta^3.
TEST(GridAdaptation, ErrorWeightConvergesToAnalytic)
{
    double d2exact = 4.0 / (3.0 * std::sqrt(3.0)) / (kDelta * kDelta);
    double d3exact = 2.0 / (kDelta * kDelta * kDelta);

    OneDimGrid coarse = makeAdaptGrid(uniformGrid(65, 0, 2e-3));
    OneDimGrid fine = makeAdaptGrid(uniformGrid(257, 0, 2e-3));
    dvector vCoarse = sampleProfile(coarse);
    dvector vFine = sampleProfile(fine);

    // p = 2: third-derivative magnitude
    double e3c = std::abs(maxWeight(coarse, vCoarse) - d3exact) / d3exact;
    double e3f = std::abs(maxWeight(fine, vFine) - d3exact) / d3exact;
    EXPECT_LT(e3f, 0.05);
    EXPECT_LT(e3f, e3c);

    // p = 1: second-derivative magnitude
    coarse.errorOrder = 1;
    fine.errorOrder = 1;
    double e2c = std::abs(maxWeight(coarse, vCoarse) - d2exact) / d2exact;
    double e2f = std::abs(maxWeight(fine, vFine) - d2exact) / d2exact;
    EXPECT_LT(e2f, 0.05);
    EXPECT_LT(e2f, e2c);
}

TEST(GridAdaptation, ErrorWeightToleratesNonuniformGrid)
{
    // Geometrically stretched grid (ratio 1.02) that still resolves the
    // front: spacing at x = kXcenter is ~2.2e-5 = delta/4.5
    dvec x(220);
    double xx = 0, h = 2e-6;
    for (int j = 0; j < 220; j++) {
        x[j] = xx;
        xx += h;
        h *= 1.02;
    }
    OneDimGrid grid = makeAdaptGrid(x);
    dvector v = sampleProfile(grid);

    double d3exact = 2.0 / (kDelta * kDelta * kDelta);
    EXPECT_NEAR(maxWeight(grid, v), d3exact, 0.4 * d3exact);
}
```

- [ ] **Step 2: Run to verify it fails to compile**

```bash
pixi run -- scons bin/unittest
```

Expected: compile error — `no member named 'computeErrorWeights' in 'OneDimGrid'` (and missing `errTol`/`errorOrder`/`errCoeff`).

- [ ] **Step 3: Add members and implementation**

In `src/grid.h`, after the `absvtol` member (line ~64), add (do not remove vtol/dvtol yet):

```cpp
    //! Local-error tolerance for grid adaptation. The estimated local
    //! representation error of each adapted component must be smaller than
    //! errTol times the component's range. Set from ConfigOptions in
    //! setOptions().
    double errTol;

    //! Spatial order p of the active convection scheme (1: firstOrderUpwind,
    //! 2: secondOrderLimited). The adaptation error estimate scales as
    //! h^(p+1) * |d^(p+1)v/dx^(p+1)|. Set in setOptions().
    int errorOrder;

    //! Leading coefficient C_p of the local error estimate
    //! E = C_p * h^(p+1) * |d^(p+1)v|. Initial values are the classical
    //! interpolation-error bounds (1/8 for p=1, 1/15 for p=2); the p=2
    //! value is calibrated so matched errTol gives matched accuracy across
    //! convection schemes. Set in setOptions().
    double errCoeff;
```

In `src/grid.h`, next to the `setOptions` declaration (line ~185), add:

```cpp
    //! Estimate the nodal magnitude of the (#errorOrder + 1)-th derivative
    //! of v, the weight in the local error estimate used by adapt().
    //! Repeated application of the nonuniform first-derivative stencils;
    //! end nodes copy the nearest interior estimate. updateValues() must
    //! have been called first.
    void computeErrorWeights(const dvector& v, dvector& W) const;
```

In `src/grid.cpp`, after `setOptions` (line ~41), add:

```cpp
void OneDimGrid::computeErrorWeights(const dvector& v, dvector& W) const
{
    W.assign(jj + 1, 0.0);
    if (jj < 2) {
        return;
    }

    // Nodal first and second derivatives from the nonuniform
    // centered-difference coefficients; end nodes copy the nearest
    // interior estimate.
    dvector d1(jj + 1), d2(jj + 1);
    for (size_t i = 1; i < jj; i++) {
        d1[i] = cfp[i] * v[i+1] + cf[i] * v[i] + cfm[i] * v[i-1];
    }
    d1[0] = d1[1];
    d1[jj] = d1[jj-1];

    for (size_t i = 1; i < jj; i++) {
        d2[i] = cfp[i] * d1[i+1] + cf[i] * d1[i] + cfm[i] * d1[i-1];
    }
    d2[0] = d2[1];
    d2[jj] = d2[jj-1];

    if (errorOrder == 1) {
        for (size_t i = 0; i <= jj; i++) {
            W[i] = std::abs(d2[i]);
        }
    } else {
        for (size_t i = 1; i < jj; i++) {
            W[i] = std::abs(cfp[i] * d2[i+1] + cf[i] * d2[i] + cfm[i] * d2[i-1]);
        }
        W[0] = W[1];
        W[jj] = W[jj-1];
    }
}
```

(Add `#include <cmath>` to grid.cpp if `std::abs` on doubles is not already available via existing includes.)

- [ ] **Step 4: Run tests to verify they pass**

```bash
pixi run -- scons bin/unittest && ./bin/unittest --gtest_filter='GridAdaptation.*'
```

Expected: both `GridAdaptation` tests PASS. If the 5% bound in `ErrorWeightConvergesToAnalytic` fails marginally, check stencil indexing before loosening anything — the repeated-stencil estimate on a uniform grid should be well within 5% at N=257 for these profiles.

- [ ] **Step 5: Run the full existing gtest suite (no regressions)**

```bash
./bin/unittest
```

Expected: all tests pass (existing `ConvectionDifferencer` etc. untouched).

- [ ] **Step 6: Commit**

```bash
git add src/grid.h src/grid.cpp test/test_gridAdaptation.cpp
git commit -m "convection: [G1] scheme-aware error-weight estimator on OneDimGrid"
```

---

### Task G2: Error-budget criterion in `adapt()` + plumbing + dampVal fix

**Files:**
- Modify: `src/grid.h` (remove `vtol_in`/`dvtol_in`/`vtol`/`dvtol` at lines 46–60; rewrite the `adapt()` doc comment at lines 150–171)
- Modify: `src/grid.cpp` (`setOptions` lines 16–41; `adapt()` insertion lines ~95–142 and removal lines ~233–300 plus the `removePoint` call site; end of `adapt()`)
- Modify: `src/readConfig.h` (lines 145–146)
- Modify: `python/ember/_ember.pxd` (line 91), `python/ember/_ember.pyx` (lines 295–296)
- Modify: `python/ember/input.py` (Grid class, line ~395)
- Test: `test/test_gridAdaptation.cpp` (new tests)

**Interfaces:**
- Consumes: `computeErrorWeights`, `errTol`/`errorOrder`/`errCoeff` from G1; existing `ConfigOptions::convectionScheme` string (`src/readConfig.h:79`).
- Produces: `ConfigOptions::errTol` (double, replaces `vtol`/`dvtol` fields); Python option `Grid.errTol` (`FloatOption(2e-3, min=0)`) forwarded as `opts.errTol`; `adapt()` behavior contract relied on by G4–G7: budget met at convergence, hysteretic removal, dampVal kept consistent.

- [ ] **Step 1: Write the failing tests**

Append to `test/test_gridAdaptation.cpp` (inside the same file; helpers go in the anonymous namespace):

```cpp
namespace {

//! Refill y from the analytic profile (mimics the solver re-solving on the
//! new grid), reset dampVal, and run one adapt() pass.
void adaptOnce(OneDimGrid& grid, std::vector<dvector>& y)
{
    y[0].resize(grid.nPoints);
    for (size_t j = 0; j < grid.nPoints; j++) {
        y[0][j] = tanhProfile(grid.x[j]);
    }
    grid.dampVal = dvec::Constant(grid.nPoints, 1e6);
    grid.updated = false;
    grid.adapt(y);
}

//! adapt() to a fixed point; returns iterations used (maxIter = no fixed point).
int adaptToConvergence(OneDimGrid& grid, std::vector<dvector>& y, int maxIter = 25)
{
    for (int i = 0; i < maxIter; i++) {
        adaptOnce(grid, y);
        if (!grid.updated) {
            return i;
        }
    }
    return maxIter;
}

//! Max-filtered interval error estimate, mirroring the criterion in adapt().
double intervalError(OneDimGrid& grid, const dvector& v, size_t j)
{
    dvector W;
    grid.computeErrorWeights(v, W);
    size_t i0 = (j == 0) ? 0 : j - 1;
    size_t i1 = std::min(j + 2, grid.jj);
    double w = 0;
    for (size_t i = i0; i <= i1; i++) {
        w = std::max(w, W[i]);
    }
    return grid.errCoeff * std::pow(grid.hh[j], grid.errorOrder + 1) * w;
}

} // namespace

TEST(GridAdaptation, MeetsErrorBudgetOnTanh)
{
    OneDimGrid grid = makeAdaptGrid(uniformGrid(17, 0, 2e-3));
    std::vector<dvector> y(1);
    grid.nAdapt = 1;

    int iters = adaptToConvergence(grid, y);
    ASSERT_LT(iters, 25) << "grid did not reach a fixed point";

    grid.updateValues();
    double vRange = mathUtils::range(y[0]);
    for (size_t j = 0; j < grid.jj; j++) {
        EXPECT_LE(intervalError(grid, y[0], j), grid.errTol * vRange * 1.000001)
            << "budget violated at interval " << j;
    }
    EXPECT_GT(grid.nPoints, 20u);
    EXPECT_LT(grid.nPoints, 500u);
}

TEST(GridAdaptation, SecondPassIsIdempotent)
{
    OneDimGrid grid = makeAdaptGrid(uniformGrid(17, 0, 2e-3));
    std::vector<dvector> y(1);
    grid.nAdapt = 1;
    ASSERT_LT(adaptToConvergence(grid, y), 25);

    adaptOnce(grid, y);
    EXPECT_FALSE(grid.updated) << "converged grid changed on an extra pass";
}

TEST(GridAdaptation, HigherOrderNeedsFewerPoints)
{
    OneDimGrid g2 = makeAdaptGrid(uniformGrid(17, 0, 2e-3));  // p=2
    OneDimGrid g1 = makeAdaptGrid(uniformGrid(17, 0, 2e-3));
    g1.errorOrder = 1;
    g1.errCoeff = 0.125;

    std::vector<dvector> y1(1), y2(1);
    g1.nAdapt = 1;
    g2.nAdapt = 1;
    ASSERT_LT(adaptToConvergence(g1, y1), 25);
    ASSERT_LT(adaptToConvergence(g2, y2), 25);
    EXPECT_GT(g1.nPoints, g2.nPoints);
}

TEST(GridAdaptation, CoarsensOverresolvedGrid)
{
    OneDimGrid grid = makeAdaptGrid(uniformGrid(513, 0, 2e-3));
    std::vector<dvector> y(1);
    grid.nAdapt = 1;

    ASSERT_LT(adaptToConvergence(grid, y), 25);
    EXPECT_LT(grid.nPoints, 300u) << "over-resolved grid was not coarsened";

    grid.updateValues();
    double vRange = mathUtils::range(y[0]);
    for (size_t j = 0; j < grid.jj; j++) {
        EXPECT_LE(intervalError(grid, y[0], j), grid.errTol * vRange * 1.000001);
    }
}

TEST(GridAdaptation, DampValRespectedWhileGridChanges)
{
    OneDimGrid grid = makeAdaptGrid(uniformGrid(17, 0, 2e-3));
    std::vector<dvector> y(1);
    grid.nAdapt = 1;

    // Damping limit hh <= dampConst*dampVal = 2e-5 in the left half only.
    // Uses a custom loop because dampVal must be re-set every pass.
    for (int i = 0; i < 25; i++) {
        y[0].resize(grid.nPoints);
        for (size_t j = 0; j < grid.nPoints; j++) {
            y[0][j] = tanhProfile(grid.x[j]);
        }
        grid.dampVal.resize(grid.nPoints);
        for (size_t j = 0; j < grid.nPoints; j++) {
            grid.dampVal[j] = (grid.x[j] < 1e-3) ? 2e-5 / grid.dampConst : 1e6;
        }
        grid.updated = false;
        grid.adapt(y);
        if (!grid.updated) {
            break;
        }
    }
    ASSERT_FALSE(grid.updated);

    grid.updateValues();
    for (size_t j = 0; j < grid.jj; j++) {
        if (grid.x[j + 1] < 1e-3) {
            EXPECT_LE(grid.hh[j], 2e-5 * 1.01)
                << "damping limit violated at j = " << j;
        }
    }
}
```

- [ ] **Step 2: Run to verify the new tests fail**

```bash
pixi run -- scons bin/unittest && ./bin/unittest --gtest_filter='GridAdaptation.*'
```

Expected: compiles (members still exist), but `MeetsErrorBudgetOnTanh`, `CoarsensOverresolvedGrid`, `DampValRespectedWhileGridChanges` FAIL — `adapt()` still uses vtol/dvtol (uninitialized in these tests may also assert/misbehave; any failure mode is acceptable evidence here as long as the tests run after Step 3).

- [ ] **Step 3: Replace the criterion in `adapt()` and wire the options**

3a. `src/readConfig.h` lines 145–146: replace

```cpp
    double vtol; //!< [grid.vtol]
    double dvtol; //!< [grid.dvtol]
```

with

```cpp
    double errTol; //!< [grid.errTol]
```

3b. `src/grid.cpp` `setOptions`: replace lines 18–19 (`vtol_in`/`dvtol_in` assignments) with

```cpp
    errTol = options.errTol;
    if (options.convectionScheme == "firstOrderUpwind") {
        errorOrder = 1;
        errCoeff = 0.125;     // 1/8: piecewise-linear representation error
    } else {
        errorOrder = 2;
        errCoeff = 1.0/15.0;  // ~quadratic representation error; calibrated in [G5]
    }
```

3c. `src/grid.h`: delete the `vtol_in`, `dvtol_in`, `vtol`, `dvtol` members (lines 46–60). Rewrite the `adapt()` doc comment (lines 150–171) to describe the new algorithm:

```cpp
    //! *Adaptation algorithm*
    //!
    //! Insertion is performed first. For each adapted component f with
    //! range(f) >= #absvtol, the local representation error of interval j
    //! is estimated as
    //!     E = C_p * hh[j]^(p+1) * max|d^(p+1) f|
    //! where p (#errorOrder) is the spatial order of the convection scheme,
    //! C_p is #errCoeff, and the derivative magnitude is taken as the
    //! maximum of the estimates (computeErrorWeights()) at the nodes within
    //! one interval of j. A point is inserted in interval j if
    //!     E > errTol * range(f)
    //! or if the damping, maximum-spacing, or uniformity criteria require
    //! one.
    //!
    //! Removal is considered next: point j is removed only if, for every
    //! component, the estimate E evaluated for the merged interval
    //! (hh[j]+hh[j-1], derivative maximum over nodes within two of j) stays
    //! below rmTol * errTol * range(f), and the damping, maximum-spacing,
    //! and uniformity criteria all permit it. Since merging roughly doubles
    //! h and E scales as h^(p+1), removal has strong built-in hysteresis.
```

3d. `src/grid.cpp` `adapt()`: at the top (after `setSize(y[0].size())`, line ~82), delete the `vtol`/`dvtol` per-component resize block (lines 84–89) and add the working `dampVal` copy:

```cpp
    // Working copy of dampVal kept consistent with the changing grid.
    // updateValues() resizes dampVal without preserving its contents, so
    // the member cannot be read safely once a point has been added or
    // removed within this pass.
    dvector dampLocal(dampVal.data(), dampVal.data() + nPoints);
```

Replace `dvector dv(jj+1);` (line ~98) with `dvector W;`.

Replace the per-variable insertion tests (lines 105–142, the `for (size_t k...)` loop body computing `dv`, `vRange`, `dvRange` and the two `if` tests) with:

```cpp
        // Consider the local error estimate for each variable v in y
        for (size_t k=0; k<nAdapt; k++) {
            dvector& v = y[k];
            double vRange = mathUtils::range(v);
            if (vRange < absvtol) {
                continue; // Ignore minor species
            }
            computeErrorWeights(v, W);

            // Local representation error estimate for interval j, using the
            // largest derivative weight within one interval of j
            size_t i0 = (j == 0) ? 0 : j - 1;
            size_t i1 = std::min(j + 2, jj);
            double w = 0;
            for (size_t i=i0; i<=i1; i++) {
                w = std::max(w, W[i]);
            }
            double E = errCoeff * pow(hh[j], errorOrder+1) * w;
            if (E > errTol*vRange) {
                insert = true;
                if (debugParameters::debugAdapt) {
                    logFile.write(format("Adapt: local error wants grid point"
                        " j = %i, k = %i; E/range = %g > %g") %
                        j % k % (E/vRange) % errTol);
                }
            }
        }
```

In the damping criterion (line ~145) and its removal twin (line ~283), replace `dampVal[j]` with `dampLocal[j]`.

In the insertion branch (line ~212, `if (insert) { ... addPoint(...); ... }`), after `addPoint(static_cast<int>(j), y);` add:

```cpp
            dampLocal.insert(dampLocal.begin() + j + 1,
                             0.5*(dampLocal[j] + dampLocal[j+1]));
```

Replace the per-variable removal tests (lines 243–280) with:

```cpp
        for (size_t k=0; k<nAdapt; k++) {
            dvector& v = y[k];
            double vRange = mathUtils::range(v);
            if (vRange < absvtol) {
                continue; // Ignore minor species
            }
            computeErrorWeights(v, W);

            // Error estimate for the interval that would result from
            // removing point j (spanning x[j-1] to x[j+1])
            size_t i0 = (j < 2) ? 0 : j - 2;
            size_t i1 = std::min(j + 2, jj);
            double w = 0;
            for (size_t i=i0; i<=i1; i++) {
                w = std::max(w, W[i]);
            }
            double E = errCoeff * pow(hh[j]+hh[j-1], errorOrder+1) * w;
            if (E > rmTol*errTol*vRange) {
                if (debugParameters::debugAdapt) {
                    logFile.write(format(
                        "Adapt: no removal - error budget. j = %i, k = %i;"
                        " E/range = %g > %g") %
                        j % k % (E/vRange) % (rmTol*errTol));
                }
                remove = false;
            }
        }
```

(The removal loop needs its own `dvector W;` if the insertion one is out of scope.)

At the `removePoint` call site in the removal loop (where `removePoint(static_cast<int>(j), y);` is called and `updated` is set, line ~335–345), add immediately after `removePoint(...)`:

```cpp
            dampLocal.erase(dampLocal.begin() + j);
```

At the end of `adapt()` (after the removal loop, before the final logging/return, line ~350), write the maintained copy back so callers see a consistent array:

```cpp
    dampVal = Eigen::Map<const dvec>(dampLocal.data(), dampLocal.size());
```

3e. `python/ember/_ember.pxd` line 91: change

```
        double vtol, dvtol, rmTol, dampConst, gridMin, gridMax
```

to

```
        double errTol, rmTol, dampConst, gridMin, gridMax
```

3f. `python/ember/_ember.pyx` lines 295–296: replace the two assignments with

```python
        opts.errTol = self.grid.errTol
```

3g. `python/ember/input.py`, in `class Grid` (before `vtol`, line ~397), add:

```python
    #: Target relative accuracy of the solution on the adapted grid. The
    #: estimated local truncation error of each state-vector component is
    #: kept below ``errTol`` times that component's range, using an error
    #: estimate that accounts for the order of the selected
    #: ``general.convectionScheme``. The same value yields similar solution
    #: accuracy under either scheme; the higher-order scheme needs fewer
    #: grid points. For high accuracy, ``errTol = 5e-4``; for minimal
    #: accuracy, ``errTol = 0.01``.
    errTol = FloatOption(2e-3, min=0)
```

(Leave `vtol`/`dvtol` defined for now — G3 handles deprecation. They are simply no longer forwarded.)

- [ ] **Step 4: Build and run the gtest suite**

```bash
pixi run -- scons bin/unittest && ./bin/unittest --gtest_filter='GridAdaptation.*'
```

Expected: all `GridAdaptation` tests PASS, including the G1 ones.

- [ ] **Step 5: Run the full test suite (C++ + Python)**

```bash
pixi run test
```

Expected: PASS. Watch for: (a) Cython rebuild picking up the pxd change; (b) `test/python/test_config.py` or `test_flame_configs.py` referencing `vtol` — if a test sets `Grid(vtol=...)`, switch it to `errTol` with a comparable resolution (`vtol=0.12` ↔ `errTol=2e-3` provisionally). Report any such edits in the commit message.

- [ ] **Step 6: Commit**

```bash
git add src/grid.h src/grid.cpp src/readConfig.h python/ember/_ember.pxd python/ember/_ember.pyx python/ember/input.py test/test_gridAdaptation.cpp
git commit -m "convection: [G2] error-budget adaptation criterion replaces vtol/dvtol; dampVal kept consistent in adapt()"
```

---

### Task G3: Deprecation UX for vtol/dvtol

**Files:**
- Modify: `python/ember/input.py` (Grid class lines ~397–405; `Config.validate()` line ~937)
- Modify: `python/ember/examples/example_laminarFlameSpeed.py` (line 26)
- Test: `test/python/test_config.py`

**Interfaces:**
- Consumes: `Grid.errTol` option (G2); `Option.isSet` (set in `input.py:178` when a user passes the kwarg).
- Produces: `Config._warnDeprecated()` — prints warnings, returns nothing; called at the top of `Config.validate()`.

- [ ] **Step 1: Write the failing tests**

Add to `test/python/test_config.py` (match the file's existing import style; `Config` and `Grid` come from `ember.input`):

```python
def test_deprecated_grid_tolerances_warn(capsys):
    conf = input.Config(input.Grid(vtol=0.1, dvtol=0.15))
    conf._warnDeprecated()
    out = capsys.readouterr().out
    assert 'deprecated' in out
    assert 'errTol' in out


def test_errtol_produces_no_warning(capsys):
    conf = input.Config(input.Grid(errTol=1e-3))
    conf._warnDeprecated()
    assert capsys.readouterr().out == ''


def test_absurd_errtol_warns(capsys):
    conf = input.Config(input.Grid(errTol=0.9))
    conf._warnDeprecated()
    assert 'errTol' in capsys.readouterr().out
```

(If the module is imported differently there — e.g. `from ember import input` vs `import ember`, follow the file. Do not add a new import style.)

- [ ] **Step 2: Run to verify failure**

```bash
pixi run -- python -m pytest test/python/test_config.py -v -k deprecated
```

Expected: FAIL — `Config` has no attribute `_warnDeprecated`. (If pytest must run via scons in this env, use `pixi run test` and read the pytest section.)

- [ ] **Step 3: Implement**

3a. `python/ember/input.py`, `class Grid`: replace the `vtol` and `dvtol` docstrings/options with:

```python
    #: Deprecated and ignored. Grid resolution is controlled by
    #: :attr:`errTol`.
    vtol = FloatOption(None, level=3)

    #: Deprecated and ignored. Grid resolution is controlled by
    #: :attr:`errTol`.
    dvtol = FloatOption(None, level=3)
```

3b. In `class Config`, add a method and call it first in `validate()` (line ~937):

```python
    def _warnDeprecated(self):
        if self.grid.vtol.isSet or self.grid.dvtol.isSet:
            print("WARNING: 'grid.vtol' and 'grid.dvtol' are deprecated and have"
                  " no effect.\n"
                  "         Grid resolution is now controlled by the local-error"
                  " tolerance 'grid.errTol';\n"
                  "         the same errTol gives similar accuracy for either"
                  " convection scheme.")
        if self.grid.errTol.value is not None and self.grid.errTol.value > 0.5:
            print("WARNING: 'grid.errTol' > 0.5 effectively disables grid"
                  " refinement.")
```

and in `validate()`:

```python
    def validate(self):
        self._warnDeprecated()
        error = False
```

3c. `python/ember/examples/example_laminarFlameSpeed.py` line 26: replace

```python
    Grid(vtol=0.1, dvtol=0.15, gridMin=5e-6, gridMax=0.001),
```

with

```python
    Grid(errTol=1.5e-3, gridMin=5e-6, gridMax=0.001),
```

(Provisional value — G6 revisits after calibration.)

- [ ] **Step 4: Run the tests**

```bash
pixi run -- python -m pytest test/python/test_config.py -v
```

Expected: all PASS, including the three new tests.

- [ ] **Step 5: Commit**

```bash
git add python/ember/input.py python/ember/examples/example_laminarFlameSpeed.py test/python/test_config.py
git commit -m "convection: [G3] deprecate vtol/dvtol with warnings; errTol docs and example update"
```

---

### Task G4: Convergence harness on the errTol ladder

**Files:**
- Modify: `test/convergence/run_convergence.py` (docstring lines ~43–49; `RUNGS` lines 95–101; `_grid_kwargs` lines ~108–113; metadata lines ~230–231)
- Create: `test/convergence/analyze_settling.py`

**Interfaces:**
- Consumes: `Grid.errTol` end-to-end (G2/G3); existing CLI (`--case`, `--scheme`, `--rung/--rungs`, `--outdir`, `--workdir`, `--retries`, `--list-rungs`).
- Produces: errTol-based `RUNGS`; result JSONs carrying `errTol` instead of `vtol`/`dvtol`; `analyze_settling.py` with `n_trajectory(case_dir) -> list[int]` (from `prof*.h5` file grid sizes) and `grid_settled(traj, tail_frac=0.25, tol_pts=3) -> bool`, used by G5/G6 analysis.

- [ ] **Step 1: Replace the rung ladder**

In `test/convergence/run_convergence.py` replace `RUNGS` (lines 95–101) with:

```python
# errTol ladder, coarse -> fine, ~2.5x per rung. gridMax co-scales as in the
# Task 2.1 vtol ladder so the far field stays bounded. Rung 5 targets
# accuracy at or beyond the old vtol rung-5 (the P2.4 failure regime).
RUNGS = [
    dict(errTol=2.0e-2, gridMax=4.0e-4),   # 0: coarsest
    dict(errTol=8.0e-3, gridMax=2.5e-4),   # 1
    dict(errTol=3.2e-3, gridMax=1.6e-4),   # 2
    dict(errTol=1.3e-3, gridMax=1.0e-4),   # 3
    dict(errTol=5.0e-4, gridMax=6.3e-5),   # 4
    dict(errTol=2.0e-4, gridMax=4.0e-5),   # 5: finest
]
```

Update `_grid_kwargs` (line ~110) to `kwargs = dict(errTol=r['errTol'], gridMax=r['gridMax'])`, the metadata dict (lines ~230–231) to record `'errTol': concrete.grid.errTol,` (drop the `dvtol` line), and the module docstring ladder description (lines ~43–49) to describe the errTol ladder. Grep the file for any other `vtol` references and update them.

- [ ] **Step 2: Create the settling analyzer**

Create `test/convergence/analyze_settling.py`:

```python
"""
Grid-settling analysis for errTol-ladder runs: extracts the grid-size
trajectory from the profNNNNNN.h5 outputs of a case work directory and
tests whether the grid reached a plateau (the P2.4 failure signature is a
monotonically growing, never-settling grid).

Usage: python analyze_settling.py <workdir> [<workdir> ...]
"""
import sys
import os
import glob
import h5py


def n_trajectory(case_dir):
    """Grid point count for each profile output, in time order."""
    files = sorted(glob.glob(os.path.join(case_dir, 'prof*.h5')))
    traj = []
    for f in files:
        with h5py.File(f, 'r') as h:
            traj.append(int(h['x'].shape[0]))
    return traj


def grid_settled(traj, tail_frac=0.25, tol_pts=3):
    """True if N varies by <= tol_pts over the last tail_frac of outputs."""
    if len(traj) < 8:
        return False
    tail = traj[int(len(traj) * (1 - tail_frac)):]
    return max(tail) - min(tail) <= tol_pts


def main():
    for case_dir in sys.argv[1:]:
        traj = n_trajectory(case_dir)
        status = 'SETTLED' if grid_settled(traj) else 'NOT SETTLED'
        print('%-60s N: %s -> %s (%d outputs)  %s' %
              (case_dir, traj[0] if traj else '-',
               traj[-1] if traj else '-', len(traj), status))


if __name__ == '__main__':
    main()
```

Verify against an existing Phase-2 run directory before trusting it (a passing old run should report SETTLED; a P2.4 failing run, if its outputs are still in `build/test/convergence-work/`, should report NOT SETTLED). If profile files are absent for old runs, verify on the smoke run from Step 3 instead.

- [ ] **Step 3: Smoke run**

```bash
pixi run -- python test/convergence/run_convergence.py --case strained --scheme secondOrderLimited --rung 2
pixi run -- python test/convergence/run_convergence.py --case strained --scheme firstOrderUpwind --rung 2
pixi run -- python test/convergence/analyze_settling.py build/test/convergence-work/strained_secondOrderLimited_rung2
```

Expected: both runs complete; result JSONs contain `errTol: 0.0032`; sol run reports SETTLED; sol N noticeably smaller than fou N at the same rung.

- [ ] **Step 4: Commit**

```bash
git add test/convergence/run_convergence.py test/convergence/analyze_settling.py
git commit -m "convection: [G4] errTol rung ladder + grid-settling analyzer"
```

---

### Task G5: Calibration of C_2 and the default errTol

**Files:**
- Modify: `src/grid.cpp` (the `errCoeff` value in `setOptions`, if calibration demands), `python/ember/input.py` (Grid.errTol default)
- Create: `test/convergence/results/calibration-notes.md` (working notes; summary goes in the spec addendum in G6)

**Interfaces:**
- Consumes: G4 harness; `plot_convergence.py` error-vs-reference computation (self-converged reference = finest completed rung per case/scheme, as in Task 2.2).
- Produces: final `errCoeff` constants and `Grid.errTol` default used by G6/G7.

- [ ] **Step 1: Run the calibration rungs (1–4, all cases, both schemes)**

```bash
for c in strained twin cylindrical; do
  for s in firstOrderUpwind secondOrderLimited; do
    pixi run -- python test/convergence/run_convergence.py --case $c --scheme $s --rungs 1 2 3 4
  done
done
```

Expected: all runs complete (rungs 1–4 are inside the previously-stable envelope). Use `--retries` default; note any retry in the notes file.

- [ ] **Step 2: Compute QoI errors and the parity ratio**

Using `plot_convergence.py` (adapt its loader if it still expects vtol keys — full plot updates happen in G6): for each case and rung, compute relative errors of `consumption_speed` and `peak_T` against the finest completed rung of the *same scheme*, then form `ratio(rung) = err_fou / err_sol` at matched errTol.

- [ ] **Step 3: Adjust C_2 if parity is off**

Decision rule: let `R` = geometric mean of `ratio` over cases and rungs (QoI = consumption_speed; sanity-check with peak_T). If `R` is within `[0.5, 2]`, keep `C_2 = 1/15`. Otherwise scale: sol error responds to `C_2` as `err_sol ∝ C_2^(2/3)`, so to multiply sol error by `R` (bringing the ratio to ~1), set

```
C_2_new = C_2 * R^(3/2)
```

Edit the `errCoeff` value in `src/grid.cpp` `setOptions` accordingly (update the comment with the calibrated value and this rationale), rebuild (`pixi run -- scons bin/unittest`, confirm `GridAdaptation.*` still pass — the unit tests do not depend on the exact constant beyond `HigherOrderNeedsFewerPoints`), and re-run rungs 2–3 for one case per scheme to confirm the ratio moved as predicted.

- [ ] **Step 4: Choose the default errTol**

From the sol error-vs-errTol curves, pick the errTol (round to one significant figure) where consumption-speed error ≈ 1e-4 and N ≈ 100 on the study cases (the owner's envelope). Set it as the `Grid.errTol` default in `python/ember/input.py` and update the docstring's "high accuracy"/"minimal accuracy" example values to bracket it (~5× tighter / ~5× looser). Record the measured (errTol, N, error) triplets in `test/convergence/results/calibration-notes.md`.

- [ ] **Step 5: Run the python config tests and commit**

```bash
pixi run -- python -m pytest test/python/test_config.py -v
git add src/grid.cpp python/ember/input.py test/convergence/results/calibration-notes.md
git commit -m "convection: [G5] calibrate errCoeff and default errTol against study cases"
```

---

### Task G6: Acceptance ladder, parity report, and the §P2.4 verdict

**Files:**
- Modify: `test/convergence/plot_convergence.py` (vtol → errTol metadata keys, axis labels)
- Modify: `docs/superpowers/specs/2026-07-05-error-based-grid-adaptation-design.md` (results addendum §A)
- Modify: `.superpowers/sdd/progress.md` (ledger entry)

**Interfaces:**
- Consumes: G4 harness + analyzer, G5 constants.
- Produces: the acceptance verdict and the errTol→QoI-error table (docs artifact).

- [ ] **Step 1: Full ladder**

```bash
for c in strained twin cylindrical; do
  for s in firstOrderUpwind secondOrderLimited; do
    pixi run -- python test/convergence/run_convergence.py --case $c --scheme $s
  done
done
pixi run -- python test/convergence/analyze_settling.py build/test/convergence-work/*rung*
```

Expected: all 36 runs complete; every run reports SETTLED.

- [ ] **Step 2: Evaluate the hard criterion**

Identify the errTol rungs at which sol accuracy meets or exceeds the old vtol rung-5 accuracy (compare QoI errors against the Task 2.2 results in `test/convergence/results/`). The hard criterion: at those rungs (expected: rungs 4–5), strained and cylindrical sol runs **complete with a settling grid**. Record per-run: completed? settled? final N; convection CVODE steps per global step vs the P2.3 table (should be bounded, not diverging).

**If any of these runs fails or fails to settle: STOP.** Write up the failure signature (grid trajectory, step counts, whether the grid settled — distinguishing spec §5's two contingencies: growth cap for grid creep vs σ-freeze for pure stiffness with a settled grid) in the spec addendum, commit what exists, and report to the owner for the contingency decision. Do not implement a contingency unprompted.

- [ ] **Step 3: Parity and documentation tables**

Update `plot_convergence.py` for the errTol metadata and produce: (a) error-vs-N per case/scheme/QoI; (b) error-vs-errTol per scheme (the parity view); (c) the errTol → measured QoI-error table (per case and scheme) destined for user docs; (d) sol/fou N-ratio at matched errTol. Parity acceptance: fou and sol QoI errors within ~2× at matched errTol across rungs 1–4 (drift per spec §2 is expected at the extremes).

- [ ] **Step 4: Write the results addendum and update the ledger**

Append `### A. Validation results (G6)` to the spec with: the hard-criterion verdict per configuration, parity table, errTol→error table, N-ratios, calibration summary from G5, and any anomalies. Update `.superpowers/sdd/progress.md` marking G1–G6 with commit ranges and the verdict.

- [ ] **Step 5: Commit**

```bash
git add test/convergence/plot_convergence.py test/convergence/results docs/superpowers/specs/2026-07-05-error-based-grid-adaptation-design.md .superpowers/sdd/progress.md
git commit -m "convection: [G6] errTol acceptance ladder — P2.4 verdict, parity + tol->error tables"
```

---

### Task G7: Baseline regression at new defaults and close-out

**Files:**
- Modify: `test/convergence/README.md` (document the errTol ladder and settling analyzer)
- Modify: `.superpowers/sdd/progress.md` (final entry)
- Possibly modify: `python/ember/input.py` / `example_laminarFlameSpeed.py` (only if regression exposes a bad default; record why)

**Interfaces:**
- Consumes: everything prior; Phase-0/1 baseline machinery (`test/convergence/run_baselines.py`, `compare_baselines.py`, `baselines-phase1/`).
- Produces: regression verdict against the amended spec §4 criterion (QoI scalars within 0.5% of phase-1 baselines; N same or fewer).

- [ ] **Step 1: Run the baseline set at new defaults**

```bash
pixi run -- python test/convergence/run_baselines.py --outdir test/convergence/results/baselines-errtol
pixi run -- python test/convergence/compare_baselines.py test/convergence/baselines-phase1 test/convergence/results/baselines-errtol
```

(Flags per each script's `--help`; run_baselines runs the standard example configs at library defaults, which now means `errTol` + `secondOrderLimited`.) Expected: all cases run to completion.

- [ ] **Step 2: Evaluate the regression criterion**

Acceptance (amended spec §4): consumption speed and peak T within 0.5% of the phase-1 baselines for every case; final N same or fewer than the phase-1 runs. `compare_baselines.py` thresholds were built for identical-config comparisons (~1e-6 floor) — read its report at the scalar level rather than pass/fail, and record the actual deltas. Known caveat from the ledger: `example_single` has ~14% termination-time scatter; judge it on scalars, not termination time, and note it.

- [ ] **Step 3: Docs and ledger close-out**

Update `test/convergence/README.md` (errTol ladder, analyzer, where the tol→error table lives). Add the final G7 ledger entry to `.superpowers/sdd/progress.md` including the regression deltas, and list the release-time items that remain (default-scheme flip in 1.7 note, `compare_baselines.py` species-mismatch loudness — unchanged from the Phase-2 pre-work list).

- [ ] **Step 4: Commit**

```bash
git add test/convergence/README.md test/convergence/results .superpowers/sdd/progress.md
git commit -m "convection: [G7] baseline regression at errTol defaults + close-out"
```
