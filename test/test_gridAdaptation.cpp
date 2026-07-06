#include "../src/grid.h"
#include "../src/mathUtils.h"
#include "../src/readConfig.h"
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

//! Pins the scheme -> (errorOrder, errCoeff) mapping that setOptions()
//! establishes (grid.cpp), so a future refactor or recalibration accidentally
//! changing these literals is caught here rather than only surfacing as a
//! convergence-study regression.
TEST(GridAdaptation, SetOptionsMapsConvectionSchemeToErrorModel)
{
    ConfigOptions opts;
    opts.errTol = 1e-4;
    opts.rmTol = 0.6;
    opts.dampConst = 7;
    opts.gridMin = 5e-7;
    opts.gridMax = 2e-4;
    opts.uniformityTol = 2.5;
    opts.absvtol = 1e-8;
    opts.centerGridMin = 1e-7;

    opts.convectionScheme = "firstOrderUpwind";
    OneDimGrid gridFOU;
    gridFOU.setOptions(opts);
    EXPECT_EQ(gridFOU.errorOrder, 1);
    EXPECT_DOUBLE_EQ(gridFOU.errCoeff, 0.125);

    opts.convectionScheme = "secondOrderLimited";
    OneDimGrid gridSOL;
    gridSOL.setOptions(opts);
    EXPECT_EQ(gridSOL.errorOrder, 2);
    EXPECT_DOUBLE_EQ(gridSOL.errCoeff, 0.0139);
}
