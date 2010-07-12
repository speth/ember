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
    ConvectionSystem();

    // The ODE function: ydot = f(t,y)
    int f(const realtype t, const sdVector& y, sdVector& ydot);
    // This uses an explicit integrator, so no Jacobian/Preconditioner is necessary

    void unroll_y(const sdVector& y); // fill in current state variables from sdvector
    void roll_y(sdVector& y) const; // fill in sdvector with current state variables
    void roll_ydot(sdVector& ydot) const; // fill in sdvector with current time derivatives

    void resize(const size_t nSpec, const size_t nPoints);

    dvector U, dUdt;
    dvector T, dTdt;
    Array2D Y, dYdt;

    double Tleft; // Temperature left boundary value
    dvector Yleft; // Mass fraction left boundary values
    double rVzero; // mass flux boundary value at j=0

    // The sum of the terms held constant for each variable in this system
    dvector Uconst;
    dvector Tconst;
    Array2D Yconst;

    // Cantera data
    CanteraGas* gas;

    // Auxiliary variables
    dvector V; // mass flux [kg/m^2*s]
    dvector rV; // (radial) mass flux (r*V) [kg/m^2*s or kg/m*rad*s]

private:
    void V2rV();
    void rV2V();

    size_t nSpec;
    size_t nVars; // == nSpec + 2

    dvector rho; // mixture density [kg/m^3]
    dvector Wmx; // mixture molecular weight [kg/kmol]
    dvector W; // species molecular weights [kg/kmol]

    // variables used internally
    dvector dUdx;
    dvector dTdx;
    Array2D dYdx;
};
