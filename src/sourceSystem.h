#pragma once

#include "mathUtils.h"
#include "sundialsUtils.h"
#include "chemistry0d.h"
#include "strainFunction.h"

class SourceSystem : public sdODE
{
    // This is the system representing the (chemical) source term at a point
public:
    SourceSystem();

    // The ODE function: ydot = f(t,y)
    int f(const realtype t, const sdVector& y, sdVector& ydot);

    // Calculate the Jacobian matrix: J = df/dy
    int Jac(const realtype t, const sdVector& y, const sdVector& ydot, sdMatrix& J);

    void unroll_y(const sdVector& y); // fill in current state variables from sdvector
    void roll_y(sdVector& y) const; // fill in sdvector with current state variables
    void roll_ydot(sdVector& ydot) const; // fill in sdvector with current time derivatives

    // Setup functions
    void resize(size_t nSpec);

    // A class that provides the strain rate and its time derivative
    StrainFunction strainFunction;

    // current state variables
    double U, dUdt; // tangential velocity
    double T, dTdt; // temperature
    dvector Y, dYdt; // species mass fractions

    // The sum of the terms held constant for each component in this system
    dvector C;

    // Cantera data
    CanteraGas* gas;

    // other parameters
    size_t nSpec;

private:


    // Physical properties
    double rho; // density [kg/m^3]
    double cp; // specific heat capacity (average) [J/kg*K]
    dvector cpSpec; // species specific heat capacity [J/mol*K]
    dvector W; // species molecular weights [kg/kmol]
    double Wmx; // mixture molecular weight [kg/mol]
    dvector hk; // species enthalpies [J/kmol]

    // Other quantities
    dvector wDot; // species net production rates [kmol/m^3*s]
    double qDot; // heat release rate per unit volume [W/m^3]
    double rhou; // density of the unburned mixture
};
