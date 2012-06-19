#pragma once

#include "mathUtils.h"

//! A Quasi-Steady-State ODE Integrator based on CHEMEQ2
class QSSIntegrator
{
public:
    QSSIntegrator();
    virtual ~QSSIntegrator() {}

    //! Initialize integrator arrays (problem size of N)
    virtual void initialize(size_t N);

    //! Set state at the start of a global timestep
    void setState(const dvec& yIn, double tStart);

    //! Take as many steps as needed to reach tf.
    //! Note that tf is relative to tstart, not absolute.
    int integrateToTime(double tf);
    int integrateOneStep(double tf); //!< Take one step toward tf without stepping past it.

    virtual void odefun(double t, const dvec& y, dvec& q, dvec& d, bool corrector=false) = 0;

    double tn; //!< Internal integrator time (relative to #tstart)
    dvec y; //!< current state vector

    bool stabilityCheck; //!< Flag to enable convergence-based stability check on timestep

    int itermax; //!< Number of corrector iterations to perform

    //! Accuracy parameter for determining the next timestep.
    double epsmin;

    //! Repeat timestep if correction is greater than #epsmax * #epsmin * #y[i] for any i.
    double epsmax;

    double dtmin; //!< Minimum timestep allowed.
    double dtmax; //!< Maximum timestep allowed.
    dvec ymin; //!< Minimum value allowed for each component
    dvec enforce_ymin; //!< Enforce the minimum value specified in #ymin for each component
    double abstol; //!< Minimum component value to consider when determining stability / step size

    int rcount; //!< Total number of timestep repeats (after failed timesteps)
    int gcount; //!< Total number of ODE function calls

    bool debug;

private:
    void getInitialStepSize(double tf); //!< Estimate the initial step size.

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
