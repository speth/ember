#pragma once

#include "mathUtils.h"

//! A Quasi-Steady-State ODE Integrator based on CHEMEQ2
class QSSIntegrator
{
public:
    QSSIntegrator(size_t Neq);

    int integrateToTime(double tf); //!< Take as many steps as needed to reach tf
    int integrateOneStep(double tf); //!< Take one step toward tf without stepping past it

    virtual void odefun(double t, const dvector& y, dvector& q, dvector& d) = 0;

    double t; //!< current integrator time
    dvector y; //!< current state vector

    bool stabilityCheck; //!< Flag to enable convergence-based stability check on timestep

    int itermax; //!< Number of corrector iterations to perform

    //! Accuracy parameter for determining the next timestep.
    double epsmin;

    //! Repeat timestep if correction is greater than #epsmax * #epsmin * #y[i] for any i.
    double epsmax;

    double dtmin; //!< Minimum timestep allowed.
    double tstart; //!< Independent variable at the start of the global timestep.

    dvector ymin; //!< minimum value allowed for each component

private:
    size_t N; //!< Number of state variables.

    dvector ym1; //!< Previous corrector iterate.
    dvector ym2; //!< Second previous corrector iterate.

    dvector scrarray; //!< temporary array.

    dvector q; //!< production rate
    dvector d; //!< loss rate
    dvector rtau; //!< ratio of timestep to timescale.
    dvector rtaus; //!< ratio of timestep to initial timescale for current timestep
    dvector y1; //!< predicted value
    dvector ys; // initial conditions for the chemical timestep
    dvector qs; // initial production rate

};
