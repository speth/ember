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

    void resize(const size_t nSpec, const size_t nPoints);

    dvector U, dUdt;
    dvector T, dTdt;
    vector<dvector> Y, dYdt;

    size_t nSpec;
    size_t nPoints;
    size_t nVars; // == nSpec + 2
    int alpha;

    double Tleft, Tright; // Temperature boundary values
    dvector Yleft, Yright; // Mass fraction boundary values

private:
    void V2rV();
    void rV2V();

    // Cantera data
    CanteraGas* gas;

    dvector rho; // mixture density [kg/m^3]
    dvector Wmx; // mixture molecular weight [kg/kmol]
    dvector rV; // (radial) mass flux (r*V) [kg/m^2*s or kg/m*rad*s]
    dvector W; // species molecular weights [kg/kmol]

    // The sum of the terms held constant for each variable in this system
    dvector Uconst;
    dvector Tconst;
    vector<dvector> Yconst;

    // variables used internally
    dvector V; // mass flux [kg/m^2*s]
    dvector dUdx;
    dvector dTdx;
    vector<dvector> dYdx;
};
