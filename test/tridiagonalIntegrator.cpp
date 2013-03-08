#include "../src/integrator.h"
#include "gtest/gtest.h"

// Solutions obtained using reference/bdf.py
const double soln[6][5] = {{0.00000000000000, 0.50000000000000, 2.00000000000000, 1.00000000000000, 0.00000000000000},
                           {0.09475912852595, 0.63130024184639, 1.61199522551290, 1.01579521276719, 0.22818950052384},
                           {0.17271644069234, 0.69321490131769, 1.34982270336185, 1.01014852549086, 0.38863616112939},
                           {0.23179316689213, 0.71176849986673, 1.17610917328091, 0.99302967646570, 0.49567553015264},
                           {0.27296869366712, 0.70706202264716, 1.05947277610782, 0.97139437597267, 0.56384834470269},
                           {0.29912090096054, 0.69144453968429, 0.97840324923541, 0.94893150017619, 0.60507643387595}};

class ConstantTridiagonalODE : public TridiagonalODE
{
public:
    ConstantTridiagonalODE() {}
    void get_A(dvec& a, dvec& b, dvec& c) {
        a = a_;
        b = b_;
        c = c_;
    }

    void get_k(dvec& k) {
        k = k_;
    }

    void resize(size_t N) {
        a_.resize(N);
        b_.resize(N);
        c_.resize(N);
        k_.resize(N);
    }

    dvec a_;
    dvec b_;
    dvec c_;
    dvec k_;
};


class TridiagonalIntegratorTest : public ::testing::Test
{
public:
    TridiagonalIntegratorTest()
        : integrator(ode)
        , y0(5)
    {
        integrator.resize(5);
        ode.a_ << 0, 1, 1, 1, 1;
        ode.b_ << -2, -2, -2, -2, -2;
        ode.c_ << 1, 1, 1, 1, 0;
        ode.k_ << 0, 0, 0, 0.2, 0.4;

        y0 << 0, 0.5, 2.0, 1.0, 0;
        integrator.set_y0(y0);
    }
protected:
    ConstantTridiagonalODE ode;
    TridiagonalIntegrator integrator;

    dvec y0;
};


TEST_F(TridiagonalIntegratorTest, Stepwise)
{
    // Check output after each timestep
    integrator.initialize(0, 0.2);
    dvec y = integrator.get_y();
    for (int j=0; j<6; j++) {
        for (int i=0; i<5; i++) {
            EXPECT_NEAR(soln[j][i], y[i], 1e-10);
        }
        integrator.step();
        y = integrator.get_y();
    }
}


TEST_F(TridiagonalIntegratorTest, FullIntegration)
{
    // Check output at final time using integrateToTime
    integrator.initialize(0, 0.2);
    integrator.integrateToTime(1.0);
    EXPECT_NEAR(integrator.get_h(), 0.2, 1e-10);
    EXPECT_NEAR(integrator.t, 1.0, 1e-10);
    dvec y = integrator.get_y();
    for (int i=0; i<5; i++) {
        EXPECT_NEAR(soln[5][i], y[i], 1e-10);
    }
}


TEST_F(TridiagonalIntegratorTest, ShortIntegration) {
    // Check for step size reduction to reach specified output time
    integrator.initialize(0, 1.3);
    integrator.integrateToTime(0.2);
    EXPECT_NEAR(integrator.t, 0.2, 1e-10);
    dvec y = integrator.get_y();
    for (int i=0; i<5; i++) {
        EXPECT_NEAR(soln[1][i], y[i], 1e-10);
    }
}
