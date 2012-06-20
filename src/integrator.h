#pragma once

#include "mathUtils.h"

class ODE
{
public:
    virtual ~ODE() {}
    // ODE function defined as ydot = f(t,y)
    virtual void f(const double t, const dvec& y, dvec& ydot) = 0;
};


class TridiagonalODE
{
public:
    TridiagonalODE() {}
    virtual ~TridiagonalODE() {}

    // ODE defined as ydot = f(t,y) = J*y + k
    virtual void get_A(dvec& a, dvec& b, dvec& c) = 0;
    virtual void get_k(dvec& k) = 0;
    virtual void resize(int N) {}
    virtual void initialize() {}
};


class Integrator
{
public:
    Integrator();
    virtual ~Integrator() {}

    // Initialization - Each of these must be called before starting integration
    virtual void set_y0(const dvec& y0);
    virtual void initialize(const double t0, const double h);

    // Accessor functions
    double get_h() const;
    double get_t() const;
    virtual const dvec& get_y() const;
    virtual const dvec& get_ydot() = 0;

    // Actually do the integration
    virtual void step() = 0; // take a single step
    virtual void integrateToTime(double tEnd) = 0;

    dvec y; // solution vector
    dvec ydot; // derivative of state vector
    double t; // current time

protected:
    double h; // timestep
    size_t N; // Dimension of y
};


class ExplicitIntegrator : public Integrator
{
    // Integrates an ODE defined as ydot = f(t,y) using Euler's method
public:
    ExplicitIntegrator(ODE& ode);
    ~ExplicitIntegrator() {}

    void set_y0(const dvec& y0);
    const dvec& get_ydot();

    // Actually do the integration
    void step(); // take a single step
    void integrateToTime(double tEnd);

private:
    ODE& myODE;
    dvec ydot; // time derivative of solution vector
};

// Integrator using 1st and 2nd order BDF for tridiagonal systems,
// with matrix inversions done using the Thomas algorithm.
class TridiagonalIntegrator : public Integrator
{
public:
    TridiagonalIntegrator(TridiagonalODE& ode);
    virtual ~TridiagonalIntegrator() {}

    // Setup
    void set_y0(const dvec& y0);
    void resize(int N_in);
    void initialize(const double t0, const double h);

    const dvec& get_ydot();

    // Actually do the integration
    void step(); // take a single step
    void integrateToTime(double tEnd);

private:
    TridiagonalODE& myODE;
    dvec a, b, c; // subdiagonal, diagonal, superdiagonal, right-hand-side
    dvec k; // constant contribution to y'
    dvec lu_b, lu_c, lu_d; // matrix elements after factorization
    dvec invDenom_; // internal
    unsigned int stepCount; // the number of steps taken since last initialization

    int N; // The system size
    dvec yprev; // previous solution
};
