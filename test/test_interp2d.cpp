#include "../src/quasi2d.h"

#include "gtest/gtest.h"


/* Data generated in Python:
    from scipy import interpolate
    x = np.array([0, 0.2, 0.8, 1.0])
    y = np.array([-2.0,-1.0,0.0,1.0,2.0])
    X, Y = np.meshgrid(x, y)
    D = (Y*Y + 3 * X) * (2 - X)
*/

class TestBilinearInterpolator : public testing::Test
{
public:
    TestBilinearInterpolator()
    {

        dmatrix data(4,5);
        data << 8.00, 2.00, 0.00, 2.00, 8.00,
                8.28, 2.88, 1.08, 2.88, 8.28,
                7.68, 4.08, 2.88, 4.08, 7.68,
                7.00, 4.00, 3.00, 4.00, 7.00;

        dvec x(4);
        x << 0.0, 0.2, 0.8, 1.0;

        dvec y(5);
        y << -2.0, -1.0, 0.0, 1.0, 2.0;

        interp.setup(data, x, y);
    }

    BilinearInterpolator interp;
};

TEST_F(TestBilinearInterpolator, AtGridPointsPlus)
{
    double eps = 1e-12;
    double tol = 1e-10;
    EXPECT_NEAR(0.00, interp.get(0.0 + eps, 0.0 + eps), tol);
    EXPECT_NEAR(2.88, interp.get(0.2 + eps, 1.0 + eps), tol);
    EXPECT_NEAR(7.68, interp.get(0.8 + eps, -2.0 + eps), tol);
    EXPECT_NEAR(2.88, interp.get(0.8 + eps, 0.0 + eps), tol);
    EXPECT_NEAR(7.00, interp.get(1.0 + eps, 2.0 + eps), tol);
}

TEST_F(TestBilinearInterpolator, AtGridPointsMinus)
{
    double eps = 1e-12;
    double tol = 1e-10;
    EXPECT_NEAR(0.00, interp.get(0.0 - eps, 0.0 - eps), tol);
    EXPECT_NEAR(2.88, interp.get(0.2 - eps, 1.0 - eps), tol);
    EXPECT_NEAR(7.68, interp.get(0.8 - eps, -2.0 - eps), tol);
    EXPECT_NEAR(2.88, interp.get(0.8 - eps, 0.0 - eps), tol);
    EXPECT_NEAR(7.00, interp.get(1.0 - eps, 2.0 - eps), tol);
}

TEST_F(TestBilinearInterpolator, IntermediatePointsX)
{
    double tol = 1e-12;
    EXPECT_NEAR(3.28, interp.get(0.4, -1.0), tol);
    EXPECT_NEAR(3.28, interp.get(0.4, 1.0), tol);
    EXPECT_NEAR(3.98, interp.get(0.75, -1.0), tol);
    EXPECT_NEAR(3.98, interp.get(0.75, 1.0), tol);
    EXPECT_NEAR(0.054, interp.get(0.01, 0.0), tol);
    EXPECT_NEAR(0.27, interp.get(0.05, 0.0), tol);
    EXPECT_NEAR(1.026, interp.get(0.19, 0.0), tol);
}

TEST_F(TestBilinearInterpolator, IntermediatePointsY)
{
    double tol = 1e-12;
    EXPECT_NEAR(5.0, interp.get(0.0, -1.5), tol);
    EXPECT_NEAR(3.5, interp.get(0.0, -1.25), tol);
    EXPECT_NEAR(1.0, interp.get(0.0, -0.5), tol);
    EXPECT_NEAR(1.0, interp.get(0.0, 0.5), tol);
    EXPECT_NEAR(4.44, interp.get(0.8, 1.1), tol);
    EXPECT_NEAR(7.32, interp.get(0.8, 1.9), tol);
}

TEST_F(TestBilinearInterpolator, IntermediatePointsXY)
{
    double tol = 1e-12;
    EXPECT_NEAR(5.29, interp.get(0.1, -1.5), tol);
    EXPECT_NEAR(2.73, interp.get(0.5, 0.5), tol);
    EXPECT_NEAR(5.69, interp.get(0.9, 1.5), tol);
}

TEST_F(TestBilinearInterpolator, Extrapolation)
{
    double tol = 1e-12;
    EXPECT_NEAR(-1.08, interp.get(-0.2, 0.0), tol);
    EXPECT_NEAR(0.24, interp.get(-0.4, 1.0), tol);
    EXPECT_NEAR(3.92, interp.get(1.2, -1.0), tol);
    EXPECT_NEAR(14.0, interp.get(0.0, 3.0), tol);
    EXPECT_NEAR(14.0, interp.get(0.0, -3.0), tol);
}
