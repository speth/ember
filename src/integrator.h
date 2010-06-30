#pragma once

#include "mathUtils.h"
#include "sundialsUtils.h"

class ODE
{
public:
    // ODE function defined as ydot = f(t,y)
    virtual void f(const double t, const dvector& y, dvector& ydot) = 0;
};


class LinearODE
{
public:
    LinearODE();
    ~LinearODE();

    // ODE defined as ydot = f(t,y) = J*y + c
    virtual void get_A(sdBandMatrix& J) = 0;
    virtual void get_C(dvector& y) = 0;
};


class Integrator
{
public:
    Integrator();
    ~Integrator();

    // Initialization - Each of these must be called before starting integration
    virtual void set_h(double dt);
    virtual void set_y0(const dvector& y0);
    virtual void set_t0(const double t0);

    // Accessor functions
    double get_h() const;
    double get_t() const;
    virtual const dvector& get_y() const;
    virtual const dvector& get_ydot() const;

    // Actually do the integration
    virtual void step() = 0; // take a single step
    virtual void step_to_time(double tEnd) = 0;

protected:
    double h; // timestep
    size_t N; // Dimension of y
    dvector y; // solution vector
    double t; // current time
};

class ExplicitIntegrator : public Integrator
{
    // Integrates an ODE defined as ydot = f(t,y) using Euler's method
public:
    ExplicitIntegrator(ODE& ode);
    ~ExplicitIntegrator();

    void set_y0(const dvector& y0);
    const dvector& get_ydot() const;

    // Actually do the integration
    void step(); // take a single step
    void step_to_time(double tEnd);

private:
    ODE& myODE;
    dvector ydot; // time derivative of solution vector
};

class BDFIntegrator : public Integrator
{
public:
    BDFIntegrator(LinearODE& ode);
    ~BDFIntegrator();

    // Setup
    void set_size(int N, int upper_bw, int lower_bw);
    void set_y0(const dvector& y0);
    void set_t0(const double t0);
    void set_dt(const double h);

    // Actually do the integration
    void step(); // take a single step
    void step_to_time(double tEnd);

private:
    LinearODE& myODE;
    sdBandMatrix* A;
    sdBandMatrix* LU;
    vector<int> p; // pivot matrix
    unsigned int stepCount; // the number of steps taken since last initialization
    int N; // The system size
    int upper_bw, lower_bw; // bandwidth of the Jacobian
    dvector yprev; // previous solution
    dvector c;
};
