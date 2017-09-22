#include "../src/diffusionSystem.h"
#include "gtest/gtest.h"

class DiffusionSystemTest : public ::testing::Test
{
public:
    DiffusionSystemTest() :
        solver(sys)
    {
        // Grid setup
        nPoints = 401;
        grid.unburnedLeft = true;
        grid.fixedBurnedVal = true;
        grid.cylindricalFlame = false;
        grid.twinFlame = false;
        grid.alpha = 0;
        grid.setSize(nPoints);
        grid.x = dvec::LinSpaced(nPoints, -2, 2);
        grid.updateValues();
        grid.updateBoundaryIndices();

        // Solver setup
        solver.resize(nPoints);

        // Diffusion system setup
        D = 0.1;
        sys.setGrid(grid);
        sys.resize(nPoints);
        sys.resetSplitConstants();
        sys.B.setConstant(1.0);
        sys.D.setConstant(D);

    }

    int nPoints;
    double D;
    OneDimGrid grid;
    DiffusionSystem sys;
    TridiagonalIntegrator solver;
};

TEST_F(DiffusionSystemTest, Constant)
{
    // Constant-valued solution should stay constant

    dvec y0 = dvec::Constant(nPoints, 5.0);
    solver.set_y0(y0);
    ASSERT_EQ(solver.y.size(), nPoints);
    solver.initialize(0.0, 0.01);
    const dvec& ydot0 = solver.get_ydot();
    for (int i = 0; i < nPoints; i++) {
        ASSERT_NEAR(ydot0[i], 0.0, 1e-11);
    }
    solver.integrateToTime(1.0);

    const dvec& ydot1 = solver.get_ydot();
    const dvec& y1 = solver.get_y();
    for (int i = 0; i < nPoints; i++) {
        ASSERT_NEAR(ydot1[i], 0.0, 1e-11);
        ASSERT_NEAR(y0[i], y1[i], 1e-11);
    }
}

TEST_F(DiffusionSystemTest, Gaussian)
{
    // A case where an analytical solution is known: An initial Gaussian pulse
    // is transformed into a broader, shorter Gaussian pulse.

    double t0 = 0.2;
    dvec y0 = exp(- grid.x * grid.x / (4 * D * t0));
    solver.set_y0(y0);
    solver.initialize(t0, 0.001);
    double t = 0.3;
    solver.integrateToTime(t);

    const dvec& y1 = solver.get_y();
    dvec y1_exact = sqrt(t0/t) * exp(-grid.x * grid.x / (4 * D * t));
    for (int i = 0; i < nPoints; i++) {
        ASSERT_NEAR(y1[i], y1_exact[i], 1e-4);
    }
}

TEST_F(DiffusionSystemTest, CylindricalCoordinates)
{
    // Exact solution in cylindrical coordinates. Same shape as the Cartesian
    // case, but decays as 1/t instead of 1/sqrt(t).

    grid.cylindricalFlame = true;
    grid.alpha = 1;
    grid.fixedBurnedVal = true;
    grid.unburnedLeft = false;

    grid.x = dvec::LinSpaced(nPoints, 0, 2);
    grid.updateValues();
    grid.updateBoundaryIndices();
    sys.setGrid(grid);

    double t0 = 0.7;
    dvec y0 = exp(- grid.x * grid.x / (4 * D * t0));
    solver.set_y0(y0);
    solver.initialize(t0, 0.001);
    double t = 0.8;
    solver.integrateToTime(t);

    const dvec& y1 = solver.get_y();
    dvec y1_exact = t0/t * exp(-grid.x * grid.x / (4 * D * t));

    for (int i = 1; i < nPoints; i++) {
        ASSERT_NEAR(y1[i], y1_exact[i], 2e-4);
    }
}

TEST_F(DiffusionSystemTest, Superposition1)
{
    // Since the diffusion equation is linear, the solution to an initial
    // condition written as the sum of two parts should be the sum of the
    // solutions to the individual parts.

    dvec y0a = cos(M_PI * grid.x);
    dvec y0b = sin(M_PI * grid.x / 4.0);

    dvec y0 = y0a + y0b;

    solver.set_y0(y0a);
    solver.initialize(0, 0.001);
    solver.integrateToTime(0.1);
    dvec y1a = solver.get_y();

    solver.set_y0(y0b);
    solver.initialize(0, 0.001);
    solver.integrateToTime(0.1);
    dvec y1b = solver.get_y();

    solver.set_y0(y0);
    solver.initialize(0, 0.001);
    solver.integrateToTime(0.1);
    dvec y1 = solver.get_y();

    for (int i = 0; i < nPoints; i++) {
        ASSERT_NEAR(y1a[i] + y1b[i], y1[i], 1e-12);
    }
}

TEST_F(DiffusionSystemTest, Superposition2)
{
    // Superposition still holds for non-uniform diffusion coefficients and
    // and prefactors.

    sys.D = 0.1 + 0.02 * grid.x;
    sys.B = 1 + 0.5 * grid.x * grid.x;

    dvec y0a = cos(M_PI * grid.x);
    dvec y0b = sin(M_PI * grid.x / 4.0);

    dvec y0 = y0a + y0b;

    solver.set_y0(y0a);
    solver.initialize(0, 0.001);
    solver.integrateToTime(0.1);
    dvec y1a = solver.get_y();

    solver.set_y0(y0b);
    solver.initialize(0, 0.001);
    solver.integrateToTime(0.1);
    dvec y1b = solver.get_y();

    solver.set_y0(y0);
    solver.initialize(0, 0.001);
    solver.integrateToTime(0.1);
    dvec y1 = solver.get_y();

    for (int i = 0; i < nPoints; i++) {
        ASSERT_NEAR(y1a[i] + y1b[i], y1[i], 1e-12);
    }
}

TEST_F(DiffusionSystemTest, Timestep)
{
    // The solution should be (nearly) independent of the timestep

    sys.D = 0.1 + 0.02 * grid.x;
    sys.B = 1 + 0.5 * grid.x * grid.x;

    double t0 = 0.2;
    dvec y0 = cos(M_PI * grid.x) + exp(- grid.x * grid.x / (4 * D * t0));

    solver.set_y0(y0);
    solver.initialize(0, 0.001);
    solver.integrateToTime(0.1);
    dvec y1a = solver.get_y();

    solver.set_y0(y0);
    solver.initialize(0, 0.0005);
    solver.integrateToTime(0.1);
    dvec y1b = solver.get_y();

    for (int i = 0; i < nPoints; i++) {
        ASSERT_NEAR(y1a[i], y1b[i], 5e-6);
    }
}
