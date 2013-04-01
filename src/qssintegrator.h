#pragma once

#include "mathUtils.h"

//! An %ODE to be solved by QssIntegrator
class QssOde
{
public:
    QssOde() {}
    virtual ~QssOde() {}
    //! Evaluate the %ODE function. The %ODE should have the form:
    //! \f[ y' = q(y) - d(y) \f]
    //! @param t Time [s] at which to evaluate the equations
    //! @param y State variables at time `t`.
    //! @param[out] q Rate at which `y` is being produced, where `q[k]` is
    //!     approximately independendent of `y[k]`.
    //! @param[out] d Rate at which `y` is being consumed / lost, where `d[k]`
    //!     is approximately linearly proportional to `y[k]`.
    //! @param corrector `true` during corrector iterations
    virtual void odefun(double t, const dvec& y, dvec& q, dvec& d, bool corrector=false) = 0;
};

//! A Quasi-Steady-State %ODE %Integrator based on CHEMEQ2
class QssIntegrator
{
public:
    QssIntegrator();
    virtual ~QssIntegrator() {}

    //! Set the ODE object to be integrated
    void setOde(QssOde* ode);

    //! Initialize integrator arrays (problem size of N)
    virtual void initialize(size_t N);

    //! Set state at the start of a global timestep
    void setState(const dvec& yIn, double tStart);

    //! Take as many steps as needed to reach `tf`.
    //! Note that `tf` is relative to #tstart, not absolute.
    int integrateToTime(double tf);

    //! Take one step toward `tf` without stepping past it.
    int integrateOneStep(double tf);

    double tn; //!< Internal integrator time (relative to #tstart)
    dvec y; //!< current state vector

    //! Flag to enable convergence-based stability check on timestep
    bool stabilityCheck;

    int itermax; //!< Number of corrector iterations to perform

    //! Accuracy parameter for determining the next timestep.
    double epsmin;

    //! Repeat timestep if correction is greater than
    //! `#epsmax * #epsmin * #y[i]` for any `i`.
    double epsmax;

    double dtmin; //!< Minimum timestep allowed.
    double dtmax; //!< Maximum timestep allowed.
    dvec ymin; //!< Minimum value allowed for each component

    //! Enforce the minimum value specified in #ymin for each component.
    dvec enforce_ymin;

    //! Minimum component value to consider when determining stability / step
    //! size.
    double abstol;

    int rcount; //!< Total number of timestep repeats (after failed timesteps)
    int gcount; //!< Total number of %ODE function calls

private:
    //! Estimate the initial step size.
    //!
    //! Strongly increasing functions(`q >> d` assumed here) use a step- size
    //! estimate proportional to the step needed for the function to reach
    //! equilibrium where as functions decreasing or in equilibrium use a
    //! stepsize estimate directly proportional to the characteristic stepsize
    //! of the function. convergence of the integration scheme is likely since
    //! the smallest estimate is chosen for the initial stepsize.
    void getInitialStepSize(double tf);

    QssOde* ode_;
    size_t N; //!< Number of state variables.

    double ts; //!< Time at the start of the current internal timestep
    double tfd; //!< Round-off parameter used to determine when integration is complete
    double dt; //!< Integrator internal timestep
    double tstart; //!< Independent variable at the start of the global timestep.

    dvec ym1; //!< Previous corrector iterate.
    dvec ym2; //!< Second previous corrector iterate.

    dvec scratch; //!< temporary array.

    dvec q; //!< production rate
    dvec d; //!< loss rate
    dvec rtau; //!< ratio of timestep to timescale.
    dvec rtaus; //!< ratio of timestep to initial timescale for current timestep
    dvec y1; //!< predicted value
    dvec ys; //!< initial conditions for the chemical timestep
    dvec qs; //!< initial production rate

    //! Flag to trigger integrator restart on the next call to
    //! #integrateToTime or #integrateOneStep.
    bool firstStep;
};
