#include "integrator.h"

using std::cout;
using std::endl;


Integrator::Integrator()
    : t(0)
{
}

void Integrator::set_h(double dt)
{
    h = dt;
}

void Integrator::set_y0(const dvector& y0)
{
    N = y0.size();
    y.assign(y0.begin(), y0.end());
}

void Integrator::set_t0(const double t0)
{
    t = t0;
}

double Integrator::get_h() const
{
    return h;
}

const dvector& Integrator::get_y() const
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

void ExplicitIntegrator::set_y0(const dvector& y0)
{
    Integrator::set_y0(y0);
    ydot.resize(N,0);
}

const dvector& ExplicitIntegrator::get_ydot()
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

void BDFIntegrator::set_y0(const dvector& y0)
{
    Integrator::set_y0(y0);
    stepCount = 0;
}

void BDFIntegrator::set_t0(const double t0)
{
    Integrator::set_t0(t0);
    stepCount = 0;
}

void BDFIntegrator::set_dt(const double h_in)
{
    Integrator::set_h(h_in);
    stepCount = 0;
}

void BDFIntegrator::initialize()
{
    stepCount = 0;
}

const dvector& BDFIntegrator::get_ydot()
{
    myODE.get_A(*A);
    sdBandMatrix& AA = *A;
    myODE.get_C(c);

    // TODO: Use something more clever for the Matrix-vector product here
    // Why doesn't sundials provide DAXPY?
    size_t N = y.size();
    ydot.resize(N);

    ydot[0] = AA(0,0)*y[0] + AA(0,1)*y[1] + c[0];
    for (size_t i=1; i<N-1; i++) {
        ydot[i] = AA(i,i-1)*y[i-1] + AA(i,i)*y[i] + AA(i,i+1)*y[i+1] + c[i];
    }
    ydot[N-1] = AA(N-1,N-2)*y[N-2] + AA(N-1,N-1)*y[N-1] + c[N-1];
    return ydot;
}

const dvector& BDFIntegrator::get_diagonal()
{
    Adiag.resize(N);
    // Assume that A has already been updated
    for (int i=0; i<N; i++) {
        Adiag[i] = (*A)(i,i);
    }
    return Adiag;
}

void BDFIntegrator::step()
{
    if (stepCount == 0) {
        yprev.assign(y.begin(), y.end());

        // Take 8 substeps using first-order BDF
        myODE.get_A(*A);
        myODE.get_C(c);
        BandCopy(A->forSundials(), LU->forSundials(), upper_bw, lower_bw);
        sdBandMatrix& M = *LU;
        for (int i=0; i<N; i++) {
            M(i,i) -= 1.0/(h/8.0);
        }
        BandGBTRF(LU->forSundials(), &p[0]);

        for (int j=0; j<8; j++) {
            // y_n -> y_n+1
            for (int i=0; i<N; i++) {
                y[i] = -y[i]/(h/8.0) - c[i];
            }
            BandGBTRS(LU->forSundials(), &p[0], &y[0]);
        }

    } else {
        if (stepCount == 1) {
            BandCopy(A->forSundials(), LU->forSundials(), upper_bw, lower_bw);
            sdBandMatrix& M = *LU;
            for (int i=0; i<N; i++) {
                M(i,i) -= 3.0/(2.0*h);
            }
            BandGBTRF(LU->forSundials(), &p[0]);
        }
        dvector tmp(y);

        for (int i=0; i<N; i++) {
            y[i] = -2*y[i]/h + yprev[i]/(2*h) - c[i];
        }
        yprev.assign(tmp.begin(), tmp.end());
        BandGBTRS(LU->forSundials(), &p[0], &y[0]);
    }

    stepCount++;
    t += h;
}

void BDFIntegrator::integrateToTime(double tEnd)
{
    // Make sure we hit tEnd after an integral number of timesteps,
    // Taking timesteps that are no longer than the given h
    int nSteps = (int) (tEnd-t)/h;
    if (t + h*nSteps != tEnd) {
        h = (tEnd - t) / (nSteps + 1);
    }
    while (t < tEnd) {
        step();
    }
}
