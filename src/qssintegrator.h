#pragma once

#include "mathUtils.h"

class QSSIntegrator
{
public:
    QSSIntegrator(size_t Neq);

    int integrateToTime(double tf); //!< Take as many steps as needed to reach tf
    int integrateOneStep(double tf); //!< Take one step toward tf without stepping past it

    void chemsp(double epsmn, double epsmx, double dtmn, double tnot,
                int itermx, int ns, dvector& ymn, double prt);

    virtual void odefun(double t, const dvector& y, dvector& q, dvector& d) = 0;
    double t; //!< current integrator time
    dvector y; // !< current state vector

    bool stabilityCheck; //!< Flag to enable convergence-based stability check on timestep

private:
    size_t N;
    double epsmin;
    double epsmax;
    double epscl;
    double dtmin;
    double tstart;
    int itermax;
    dvector ymin;
};
