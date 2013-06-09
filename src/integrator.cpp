#include "integrator.h"

#include "debugUtils.h"

Integrator::Integrator()
    : t(0)
{
}

void Integrator::set_y0(const dvec& y0)
{
    N = y0.rows();
    y = y0;
}

void Integrator::initialize(const double t0, const double h_)
{
    t = t0;
    h = h_;
}

double Integrator::get_h() const
{
    return h;
}

const dvec& Integrator::get_y() const
{
    return y;
}

double Integrator::get_t() const
{
    return t;
}

void Integrator::integrateToTime(double tEnd)
{
    // Make sure we hit tEnd after an integral number of timesteps,
    // Taking timesteps that are no longer than the given h
    int nSteps = static_cast<int>((tEnd-t)/h + 1e-5);
    nSteps = std::max(nSteps, 1);
    h = (tEnd - t) / nSteps;
    for (int i = 0; i < nSteps; i++) {
        step();
    }
}

ExplicitIntegrator::ExplicitIntegrator(ODE& ode)
    : myODE(ode)
{
}

void ExplicitIntegrator::set_y0(const dvec& y0)
{
    Integrator::set_y0(y0);
    ydot.resize(N);
    ydot.setZero();
}

const dvec& ExplicitIntegrator::get_ydot()
{
    return ydot;
}

void ExplicitIntegrator::step()
{
    // Integrate one time step using the explicit Euler's method
    myODE.f(t, y, ydot);
    y += h * ydot;
    t += h;
}

// *******************************
// * class TridiagonalIntegrator *
// *******************************

TridiagonalIntegrator::TridiagonalIntegrator(TridiagonalODE& ode)
    : myODE(ode)
    , stepCount(0)
{
}


void TridiagonalIntegrator::resize(size_t N_in)
{
    N = static_cast<int>(N_in);

    y.resize(N);
    yprev.resize(N);

    myODE.resize(N);

    a.setZero(N);
    b.setZero(N);
    c.setZero(N);

    lu_b.resize(N);
    lu_c.resize(N);
    lu_d.resize(N);
    invDenom_.resize(N);
}


void TridiagonalIntegrator::set_y0(const dvec& y0)
{
    Integrator::set_y0(y0);
    stepCount = 0;
}

void TridiagonalIntegrator::initialize(const double t0, const double h_)
{
    Integrator::initialize(t0, h_);
    stepCount = 0;
    myODE.initialize();
}

const dvec& TridiagonalIntegrator::get_ydot()
{
    myODE.get_A(a, b, c);
    myODE.get_k(k);

    ydot.resize(N);
    ydot[0] = b[0]*y[0] + c[0]*y[1] + k[0];
    for (int i=1; i<N-1; i++) {
        ydot[i] = a[i]*y[i-1] + b[i]*y[i] + c[i]*y[i+1] + k[i];
    }
    ydot[N-1] = a[N-1]*y[N-2] + b[N-1]*y[N-1] + k[N-1];

    return ydot;
}

void TridiagonalIntegrator::step()
{
    if (stepCount == 0) {
        yprev = y; // current value of y becomes y_(n-1)

        // Get ODE coefficients
        myODE.get_A(a, b, c);
        myODE.get_k(k);

        // Modify diagonal elements according to 1st order BDF
        lu_b = b - 1.0/(h/8.0);

        // Compute Thomas coefficients
        lu_c[0] = c[0] / lu_b[0];
        for (int i=1; i<N; i++) {
            invDenom_[i] = 1.0 / (lu_b[i] - lu_c[i-1] * a[i]);
            lu_c[i] = c[i] * invDenom_[i];
        }

        // Take 8 substeps using first-order BDF
        for (int j=0; j<8; j++) {
            // RHS for 1st order BDF
            y = -y/(h/8.0) - k;

            // Thomas coefficient depending on y
            lu_d[0] = y[0] / lu_b[0];
            for (int i=1; i<N; i++) {
                lu_d[i] = (y[i] - lu_d[i-1]*a[i]) * invDenom_[i];
            }

            // Back substitution
            y[N-1] = lu_d[N-1];
            for (int i=N-2; i>=0; i--) {
                y[i] = lu_d[i] - lu_c[i] * y[i+1];
            }

            assert(mathUtils::notnan(y));
        }

    } else {
        if (stepCount == 1) {
            // Modify diagonal elements according to 2nd order BDF
            lu_b = b - 3.0/(2.0*h);

            // Compute Thomas coefficients
            lu_c[0] = c[0] / lu_b[0];
            for (int i=1; i<N; i++) {
                invDenom_[i] = 1.0 / (lu_b[i] - lu_c[i-1] * a[i]);
                lu_c[i] = c[i] * invDenom_[i];
            }
        }
        dvec ynm1 = y; // Current value of y becomes y_(n-1)
        y = -2 * ynm1 / h + yprev / (2*h) - k; // RHS for 2nd order BDF
        yprev = ynm1; // y_(n-1) will be y_(n-2) next time

        // Compute Thomas coefficients
        lu_d[0] = y[0] / lu_b[0];
        for (int i=1; i<N; i++) {
            lu_d[i] = (y[i] - lu_d[i-1]*a[i]) * invDenom_[i];
        }

        // Back substitution
        y[N-1] = lu_d[N-1];
        for (int i=N-2; i>=0; i--) {
            y[i] = lu_d[i] - lu_c[i] * y[i+1];
        }
        assert(mathUtils::notnan(y));
    }

    stepCount++;
    t += h;
}

