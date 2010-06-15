#pragma once

#include "mathUtils.h"

class ODE
{
public:
    // ODE function defined as ydot = f(t,y)
    virtual void f(const double t, const dvector& y, dvector& ydot) = 0;
};

class ExplicitIntegrator
{
    // Integrates an ODE defined as ydot = f(t,y) using Euler's method
public:
    ExplicitIntegrator(ODE& ode);
    ~ExplicitIntegrator();

    // Initialization - Each of these must be called before starting integration
    void set_timestep(double dt);
    void set_inital_condition(const dvector& y0);

    // Accessor functions
    double get_timestep() const;
    const dvector& get_y() const;
    const dvector& get_ydot() const;

    // Actually do the integration
    void step(); // take a single step
    void step_to_time(double t_end);

private:
    ODE& myODE;
    double timestep;

    size_t N; // Dimension of y
    dvector y; // solution vector
    dvector ydot; // time derivative of solution vector
    dvector t; // current time
};
