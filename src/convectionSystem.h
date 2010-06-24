#pragma once

#include "mathUtils.h"
#include "sundialsUtils.h"
#include "strainFunction.h"
#include "grid.h"
#include "chemistry0d.h"


class ConvectionSystem  : public sdODE, public GridBased
{
    // This is the system representing convection of all state variables in the domain.
public:
    // The ODE function: ydot = f(t,y)
    int f(const realtype t, const sdVector& y, sdVector& ydot);
    // This uses an explicit integrator, so no Jacobian/Preconditioner is necessary

    void unroll_y(const sdVector& y); // fill in current state variables from sdvector
    void roll_y(sdVector& y) const; // fill in sdvector with current state variables
    void roll_ydot(sdVector& ydot) const; // fill in sdvector with current time derivatives

    dvector U, dUdt;
    dvector T, dTdt;
    vector<dvector> Y, dYdt;

    size_t nSpec;
    size_t nPoints;
    size_t nVars; // == nSpec + 2

private:
    // Cantera data
    CanteraGas* gas;

    dvector drhodt;
    dvector Wmx; // mixture molecular weight [kg/kmol]
    dvector rV; // (radial) mass flux (r*V) [kg/m^2*s or kg/m*rad*s]
    dvector W; // species molecular weights [kg/kmol]

    // The sum of the terms held constant for this system
    dvector constantTerm;
};
