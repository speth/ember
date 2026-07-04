#include "../src/convectionDifferencer.h"
#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <random>
#include <vector>

using Scheme = ConvectionDifferencer::Scheme;

namespace {

//! Build a minimal planar OneDimGrid from a set of node coordinates.
OneDimGrid makeGrid(const dvec& x)
{
    OneDimGrid grid;
    grid.alpha = 0;
    grid.beta = 1;
    grid.setSize(x.size());
    grid.x = x;
    grid.updateValues();
    return grid;
}

//! Independent copy of the legacy 1st-order upwind convection loop, used as the
//! bit-identical reference. Forward difference where the velocity is negative
//! or at the left boundary node, backward difference otherwise.
dvec legacyUpwind(const dvec& y, const dvec& v, const OneDimGrid& grid)
{
    const dvec& hh = grid.hh;
    dvec d = dvec::Zero(grid.nPoints);
    for (size_t j = 0; j < grid.jj; j++) {
        if (v[j] < 0 || j == 0) {
            d[j] = (y[j+1] - y[j]) / hh[j];
        } else {
            d[j] = (y[j] - y[j-1]) / hh[j-1];
        }
    }
    return d;
}

//! Node coordinates for a smoothly, strongly nonuniform grid on [-L, L]. The
//! spacing density is exp(sum of a few sine modes), so adjacent spacings vary
//! smoothly (representative of the production adaptive grid) while reaching
//! adjacent-spacing ratios of order `uniformityTol` (~2) at coarse resolution.
dvec makeStretchedX(size_t N, double L,
                    const std::array<double, 3>& amps,
                    const std::array<double, 3>& phases)
{
    dvec dens(N);
    for (size_t i = 0; i < N; i++) {
        double xi = double(i) / (N - 1);
        double e = 0.0;
        for (int k = 0; k < 3; k++) {
            e += amps[k] * std::sin(2 * M_PI * (k + 1) * xi + phases[k]);
        }
        dens[i] = std::exp(e);
    }
    dvec x(N);
    x[0] = 0.0;
    for (size_t i = 1; i < N; i++) {
        x[i] = x[i-1] + 0.5 * (dens[i-1] + dens[i]); // trapezoidal integral
    }
    double span = x[N-1];
    for (size_t i = 0; i < N; i++) {
        x[i] = x[i] / span * (2 * L) - L;
    }
    return x;
}

} // anonymous namespace

// Constant field -> derivative is identically zero for both schemes, for a
// mixed-sign velocity field.
TEST(ConvectionDifferencer, ConstantField)
{
    dvec x(6);
    x << 0.0, 0.3, 0.5, 1.1, 1.4, 2.0; // nonuniform
    OneDimGrid grid = makeGrid(x);

    dvec y = dvec::Constant(6, 3.7);
    dvec v(6);
    v << 1.0, -1.0, 1.0, -1.0, 1.0, -1.0; // mixed sign

    ConvectionDifferencer diff;
    dvec d = dvec::Zero(6);
    for (Scheme s : {Scheme::FirstOrderUpwind, Scheme::SecondOrderLimited}) {
        diff.setScheme(s);
        diff.resize(6);
        diff.computeDerivatives(y, v, grid, d);
        for (size_t j = 0; j < grid.jj; j++) {
            EXPECT_NEAR(d[j], 0.0, 1e-13) << "j = " << j;
        }
    }
}

// Linear field -> exact derivative for both schemes, for a mixed-sign velocity.
// FirstOrderUpwind is exact at nodes 0..jj-1; SecondOrderLimited is exact at
// interior nodes 1..jj-1 and sets node 0 to zero.
TEST(ConvectionDifferencer, LinearField)
{
    dvec x(7);
    x << 0.0, 0.2, 0.5, 0.9, 1.0, 1.6, 2.1; // nonuniform
    OneDimGrid grid = makeGrid(x);

    double a = -1.3, b = 0.4;
    dvec y = a * x + b;

    // Mixed-sign interior velocity. The outermost node slopes are pinned to
    // zero (sigma[0] = sigma[jj] = 0), so a reconstruction that reaches across
    // one of those faces is only locally 1st-order by design (spec 3.2). That
    // happens for node 1 when v > 0 and for node jj-1 when v < 0; the velocity
    // here keeps those two nodes on their interior-facing branch so every
    // interior node is genuinely 2nd-order and therefore linear-exact.
    dvec v(7);
    v << 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0;

    ConvectionDifferencer diff;
    dvec d = dvec::Zero(7);

    // FirstOrderUpwind is linear-exact at nodes 0..jj-1 for any velocity.
    diff.setScheme(Scheme::FirstOrderUpwind);
    diff.resize(7);
    diff.computeDerivatives(y, v, grid, d);
    for (size_t j = 0; j < grid.jj; j++) {
        EXPECT_NEAR(d[j], a, 1e-10) << "FirstOrderUpwind, j = " << j;
    }

    diff.setScheme(Scheme::SecondOrderLimited);
    diff.computeDerivatives(y, v, grid, d);
    EXPECT_NEAR(d[0], 0.0, 1e-13) << "SecondOrderLimited node 0 must be 0";
    for (size_t j = 1; j < grid.jj; j++) {
        EXPECT_NEAR(d[j], a, 1e-10) << "SecondOrderLimited, j = " << j;
    }
}

// FirstOrderUpwind must reproduce the legacy loop bit-for-bit for random data
// and random velocities on a nonuniform grid.
TEST(ConvectionDifferencer, FirstOrderMatchesLegacy)
{
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> ud(-5.0, 5.0);
    std::uniform_real_distribution<double> uh(0.1, 1.0);

    const size_t N = 30;
    dvec x(N);
    x[0] = 0.0;
    for (size_t i = 1; i < N; i++) {
        x[i] = x[i-1] + uh(rng); // random nonuniform spacing
    }
    OneDimGrid grid = makeGrid(x);

    ConvectionDifferencer diff;
    diff.setScheme(Scheme::FirstOrderUpwind);
    diff.resize(N);

    for (int trial = 0; trial < 20; trial++) {
        dvec y(N), v(N);
        for (size_t i = 0; i < N; i++) {
            y[i] = ud(rng);
            v[i] = ud(rng);
        }
        dvec d = dvec::Zero(N);
        diff.computeDerivatives(y, v, grid, d);
        dvec ref = legacyUpwind(y, v, grid);
        for (size_t j = 0; j < grid.jj; j++) {
            EXPECT_EQ(d[j], ref[j]) << "trial " << trial << ", j = " << j;
        }
    }
}

// SecondOrderLimited achieves ~2nd-order convergence in the interior L2 norm on
// smoothly, strongly nonuniform grids under refinement, for velocity fields
// with sign changes. Randomized but seeded (deterministic).
TEST(ConvectionDifferencer, SecondOrderConvergence)
{
    std::mt19937 rng(20260704);
    std::uniform_real_distribution<double> uphase(0.0, 2 * M_PI);
    std::uniform_real_distribution<double> uamp(0.7, 1.3);
    const std::array<double, 3> baseAmp = {2.2, 1.1, 0.6};

    const double w = 0.5, L = 2.5;
    const std::vector<size_t> Ns = {81, 161, 321, 641};

    for (int c = 0; c < 4; c++) {
        std::array<double, 3> phases = {uphase(rng), uphase(rng), uphase(rng)};
        std::array<double, 3> amps = {baseAmp[0] * uamp(rng),
                                      baseAmp[1] * uamp(rng),
                                      baseAmp[2] * uamp(rng)};

        std::vector<double> errs;
        double maxRatio = 1.0;
        for (size_t N : Ns) {
            dvec x = makeStretchedX(N, L, amps, phases);
            OneDimGrid grid = makeGrid(x);

            dvec y = (x / w).tanh();
            dvec cosh = (x / w).cosh();
            dvec exact = (1.0 / w) * (cosh * cosh).inverse(); // sech^2(x/w)/w
            dvec v = (2.7 * x + double(c)).sin();             // sign changes

            ConvectionDifferencer diff;
            diff.setScheme(Scheme::SecondOrderLimited);
            diff.resize(N);
            dvec d = dvec::Zero(N);
            diff.computeDerivatives(y, v, grid, d);

            double se = 0.0;
            size_t cnt = 0;
            for (size_t j = 1; j < grid.jj; j++) {
                se += (d[j] - exact[j]) * (d[j] - exact[j]);
                cnt++;
                double r = grid.hh[j] / grid.hh[j-1];
                maxRatio = std::max(maxRatio, std::max(r, 1.0 / r));
            }
            errs.push_back(std::sqrt(se / cnt));
        }

        // Observed order between the two finest grids (asymptotic regime).
        size_t n = Ns.size();
        double order = std::log(errs[n-2] / errs[n-1])
                     / std::log(double(Ns[n-1]) / double(Ns[n-2]));
        EXPECT_GE(order, 1.9) << "case " << c << ": observed order " << order
                              << ", max adjacent spacing ratio " << maxRatio;
    }
}

// No-new-extrema (boundedness): on monotone data the limited reconstruction
// must not introduce over/undershoots, so the derivative keeps a consistent
// sign. Tested on a strongly nonuniform grid with a steep monotone front.
TEST(ConvectionDifferencer, NoNewExtrema)
{
    std::mt19937 rng(777);
    std::uniform_real_distribution<double> uu(-1.0, 1.0);

    const size_t N = 41;
    dvec x(N);
    x[0] = 0.0;
    double h = 1.0;
    for (size_t i = 0; i < N - 1; i++) {
        if (i) {
            h *= std::pow(2.0, uu(rng)); // adjacent ratio in [1/2, 2]
        }
        x[i+1] = x[i] + h;
    }
    x = x / x[N-1] * 4.0 - 2.0;
    OneDimGrid grid = makeGrid(x);

    dvec yInc = (x / 0.05).tanh(); // monotone increasing, steep front

    ConvectionDifferencer diff;
    diff.setScheme(Scheme::SecondOrderLimited);
    diff.resize(N);
    dvec d = dvec::Zero(N);

    for (double vs : {1.0, -1.0}) {
        dvec v = dvec::Constant(N, vs);
        diff.computeDerivatives(yInc, v, grid, d);
        for (size_t j = 1; j < grid.jj; j++) {
            EXPECT_GE(d[j], -1e-12)
                << "monotone-increasing data gave a negative slope at j = " << j
                << " (v = " << vs << ")";
        }
    }

    dvec yDec = -yInc; // monotone decreasing
    for (double vs : {1.0, -1.0}) {
        dvec v = dvec::Constant(N, vs);
        diff.computeDerivatives(yDec, v, grid, d);
        for (size_t j = 1; j < grid.jj; j++) {
            EXPECT_LE(d[j], 1e-12)
                << "monotone-decreasing data gave a positive slope at j = " << j
                << " (v = " << vs << ")";
        }
    }
}

// Upwind-branch selection (FirstOrderUpwind): a mixed-sign velocity field must
// pick the correct one-sided difference at each node. Hand-computed small case.
TEST(ConvectionDifferencer, UpwindBranchSelection)
{
    dvec x(5);
    x << 0.0, 1.0, 3.0, 6.0, 10.0; // hh = 1, 2, 3, 4
    OneDimGrid grid = makeGrid(x);

    dvec y(5);
    y << 0.0, 1.0, 4.0, 9.0, 16.0;
    dvec v(5);
    v << 1.0, -2.0, 3.0, -1.0, 5.0;

    ConvectionDifferencer diff;
    diff.setScheme(Scheme::FirstOrderUpwind);
    diff.resize(5);
    dvec d = dvec::Zero(5);
    diff.computeDerivatives(y, v, grid, d);

    // j=0: forced forward  -> (y1-y0)/hh0 = 1/1
    // j=1: v<0  -> forward  -> (y2-y1)/hh1 = 3/2
    // j=2: v>0  -> backward -> (y2-y1)/hh1 = 3/2
    // j=3: v<0  -> forward  -> (y4-y3)/hh3 = 7/4
    EXPECT_DOUBLE_EQ(d[0], 1.0);
    EXPECT_DOUBLE_EQ(d[1], 1.5);
    EXPECT_DOUBLE_EQ(d[2], 1.5);
    EXPECT_DOUBLE_EQ(d[3], 1.75);
}

// Upwind-branch selection (SecondOrderLimited): the reconstruction must use the
// upwind side per sign(v[j]), with the correct spacing in each face offset.
// Verified against the explicit van Albada reconstruction at an interior node.
TEST(ConvectionDifferencer, SecondOrderUpwindBranch)
{
    dvec x(5);
    x << 0.0, 1.0, 3.0, 6.0, 10.0; // hh = 1, 2, 3, 4
    OneDimGrid grid = makeGrid(x);

    dvec y(5);
    y << 0.0, 2.0, 3.0, 7.0, 8.0;

    auto va = [](double sm, double sp) {
        return (sm * sp > 0) ? sm * sp * (sm + sp) / (sm * sm + sp * sp) : 0.0;
    };
    const double hh1 = 2.0, hh2 = 3.0, hh3 = 4.0;
    double sig1 = va((y[1]-y[0]) / 1.0, (y[2]-y[1]) / hh1);
    double sig2 = va((y[2]-y[1]) / hh1, (y[3]-y[2]) / hh2);
    double sig3 = va((y[3]-y[2]) / hh2, (y[4]-y[3]) / hh3);
    double dlj2 = 0.5 * (x[3] - x[1]);

    ConvectionDifferencer diff;
    diff.setScheme(Scheme::SecondOrderLimited);
    diff.resize(5);
    dvec d = dvec::Zero(5);

    // v > 0 at node 2: reconstruct from the left/own side.
    dvec vpos = dvec::Constant(5, 1.0);
    diff.computeDerivatives(y, vpos, grid, d);
    double posExpected = ((y[2] + sig2 * hh2 / 2) - (y[1] + sig1 * hh1 / 2)) / dlj2;
    EXPECT_NEAR(d[2], posExpected, 1e-12);

    // v < 0 at node 2: reconstruct from the right/own side.
    dvec vneg = dvec::Constant(5, -1.0);
    diff.computeDerivatives(y, vneg, grid, d);
    double negExpected = ((y[3] - sig3 * hh2 / 2) - (y[2] - sig2 * hh1 / 2)) / dlj2;
    EXPECT_NEAR(d[2], negExpected, 1e-12);

    // The two branches must genuinely differ, so branch selection matters.
    EXPECT_NE(posExpected, negExpected);
}

#ifndef NDEBUG
// The SecondOrderLimited kernel requires at least 3 grid points. (Only active
// in debug builds where asserts are enabled.)
TEST(ConvectionDifferencer, AssertsMinimumSize)
{
    dvec x(2);
    x << 0.0, 1.0;
    OneDimGrid grid = makeGrid(x);
    dvec y(2);
    y << 0.0, 1.0;
    dvec v = dvec::Constant(2, 1.0);
    dvec d = dvec::Zero(2);

    ConvectionDifferencer diff;
    diff.setScheme(Scheme::SecondOrderLimited);
    diff.resize(2);
    EXPECT_DEATH(diff.computeDerivatives(y, v, grid, d), "");
}
#endif
