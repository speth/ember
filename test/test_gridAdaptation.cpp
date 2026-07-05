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
