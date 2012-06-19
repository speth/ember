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

void ExplicitIntegrator::integrateToTime(double tEnd)
{
    // Make sure we hit tEnd after an integral number of timesteps
    int nSteps = (int) (tEnd-t)/h;
    if (t + h*nSteps != tEnd) {
        h = (tEnd - t) / (nSteps + 1);
    }
    while (t < tEnd) {
        step();
    }
}

// ***********************
// * class BDFIntegrator *
// ***********************

BDFIntegrator::BDFIntegrator(LinearODE& ode)
    : myODE(ode)
    , A(NULL)
    , LU(NULL)
    , stepCount(0)
{
}

BDFIntegrator::~BDFIntegrator()
{
    delete LU;
    delete A;
}

void BDFIntegrator::resize(int N_in, int upper_bw_in, int lower_bw_in)
{
    N = N_in;
    upper_bw = upper_bw_in;
    lower_bw = lower_bw_in;

    delete LU;
    delete A;
    LU = new sdBandMatrix(N, upper_bw, lower_bw);
    A = new sdBandMatrix(N, upper_bw, lower_bw);
    p.resize(N);
    y.resize(N);

    myODE.resize(N);
}

void BDFIntegrator::set_y0(const dvec& y0)
{
    Integrator::set_y0(y0);
    stepCount = 0;
}

void BDFIntegrator::initialize(const double t0, const double h_)
{
    Integrator::initialize(t0, h_);
    stepCount = 0;
    myODE.initialize();
}

const dvec& BDFIntegrator::get_ydot()
{
    myODE.get_A(*A);
    sdBandMatrix& AA = *A;
    myODE.get_C(c);

    // TODO: Use something more clever for the Matrix-vector product here
    // Why doesn't sundials provide DAXPY?
    size_t N = y.rows();
    ydot.resize(N);

    ydot[0] = AA(0,0)*y[0] + AA(0,1)*y[1] + c[0];
    for (size_t i=1; i<N-1; i++) {
        ydot[i] = AA(i,i-1)*y[i-1] + AA(i,i)*y[i] + AA(i,i+1)*y[i+1] + c[i];
    }
    ydot[N-1] = AA(N-1,N-2)*y[N-2] + AA(N-1,N-1)*y[N-1] + c[N-1];
    return ydot;
}

void BDFIntegrator::step()
{
    if (stepCount == 0) {
        yprev = y;

        // Take 8 substeps using first-order BDF
        myODE.get_A(*A);
        myODE.get_C(c);
        BandCopy(A->forSundials(), LU->forSundials(), upper_bw, lower_bw);
        sdBandMatrix& M = *LU;
        for (int i=0; i<N; i++) {
            M(i,i) -= 1.0/(h/8.0);
        }

        int ierr = BandGBTRF(LU->forSundials(), &p[0]);
        assert(ierr == 0);
        assert(mathUtils::notnan(y));

        for (int j=0; j<8; j++) {
            // y_n -> y_n+1
            y = -y / (h/8.0) - c;
            BandGBTRS(LU->forSundials(), &p[0], y.data());
            assert(mathUtils::notnan(y));
        }

    } else {
        if (stepCount == 1) {
            BandCopy(A->forSundials(), LU->forSundials(), upper_bw, lower_bw);
            sdBandMatrix& M = *LU;
            for (int i=0; i<N; i++) {
                M(i,i) -= 3.0/(2.0*h);
            }
            int ierr = BandGBTRF(LU->forSundials(), &p[0]);
            assert(ierr == 0);
        }
        dvec tmp(y);
        y = -2*y/h + yprev/(2*h) - c;
        yprev = tmp;
        BandGBTRS(LU->forSundials(), &p[0], y.data());
        assert(mathUtils::notnan(y));
    }

    stepCount++;
    t += h;
}

void BDFIntegrator::integrateToTime(double tEnd)
{
    // Make sure we hit tEnd after an integral number of timesteps,
    // Taking timesteps that are no longer than the given h
    int nSteps = static_cast<int>((tEnd-t)/h);
    if (t + h*nSteps != tEnd) {
        h = (tEnd - t) / (nSteps + 1);
    }
    while (t < tEnd) {
        step();
    }
}


// *******************************
// * class TridiagonalIntegrator *
// *******************************

TridiagonalIntegrator::TridiagonalIntegrator(TridiagonalODE& ode)
    : myODE(ode)
    , stepCount(0)
{
}


void TridiagonalIntegrator::resize(int N_in)
{
    N = N_in;

    y.resize(N);
    y_in.resize(N);
    yprev.resize(N);

    myODE.resize(N);

    a.resize(N);
    b.resize(N);
    c.resize(N);

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

void TridiagonalIntegrator::integrateToTime(double tEnd)
{
    // Make sure we hit tEnd after an integral number of timesteps,
    // Taking timesteps that are no longer than the given h
    int nSteps = static_cast<int>((tEnd-t)/h);
    if (t + h*nSteps != tEnd) {
        h = (tEnd - t) / (nSteps + 1);
    }
    while (t < tEnd) {
        step();
    }
}
