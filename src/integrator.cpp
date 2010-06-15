#include "integrator.h"

ExplicitIntegrator::ExplicitIntegrator(ODE& ode)
{
    myODE = ode;
    t = 0;
}

void ExplicitIntegrator::set_timestep(double dt)
{
    timestep = dt;
}

void ExplicitIntegrator::set_inital_condition(const dvector& y0)
{
    N = y0.size();

    y.resize(N);
    ydot.resize(N,0);

    std::copy(y0.begin(), y0.end(), y.begin());
}

double ExplicitIntegrator::get_timestep() const
{
    return timestep;
}

const dvector& get_y() const
{
    return y;
}

const dvector& get_ydot() const
{
    return ydot;
}

void ExplicitIntegrator::step()
{
    // Integrate one time step using the explicit Euler's method
    myODE.f(t, y, ydot);
    y += timestep * ydot;
    t += timestep;
}

void ExplicitIntegrator::step_to_time(double t_end)
{
    while (t < tEnd) {
        step();
    }
}
