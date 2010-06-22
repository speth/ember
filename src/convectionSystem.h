#pragma once

#include "mathUtils.h"
#include "sundialsUtils.h"
#include "strainFunction.h"

class ConvectionSystem : public sdODE, public GridBased
{
    // This is the system representing convection of all state variables in the domain.
public:
    // The ODE function: ydot = f(t,y)
    int f(realtype t, sdVector& y, sdVector& ydot);
    // This uses an explicit integrator, so no Jacobian/Preconditioner is necessary

private:
    dvector drhodt;
    dvector Wmx; // mixture molecular weight [kg/kmol]
    dvector rV; // (radial) mass flux (r*V) [kg/m^2*s or kg/m*rad*s]
    dvector W; // species molecular weights [kg/kmol]

    // The sum of the terms held constant for this system
    dvector constantTerm;
};
