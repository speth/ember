#pragma once

#include "mathUtils.h"
#include "sundialsUtils.h"
#include "chemistry0d.h"
#include "strainFunction.h"
#include "../adapchem/ckcompat.h"
#include "perfTimer.h"

#include <boost/shared_ptr.hpp>

class SourceSystem : public sdODE
{
    // This is the system representing the (chemical) source term at a point
public:
    SourceSystem();

    // The ODE function: ydot = f(t,y)
    int f(const realtype t, const sdVector& y, sdVector& ydot);

    // Calculate the Jacobian matrix: J = df/dy
    int denseJacobian(const realtype t, const sdVector& y, const sdVector& ydot, sdMatrix& J);

    // A simpler finite difference based Jacobian
    int fdJacobian(const realtype t, const sdVector& y, const sdVector& ydot, sdMatrix& J);

    void unroll_y(const sdVector& y); // fill in current state variables from sdvector
    void roll_y(sdVector& y) const; // fill in sdvector with current state variables
    void roll_ydot(sdVector& ydot) const; // fill in sdvector with current time derivatives

    // Setup functions
    void resize(size_t nSpec);
    void initialize();

    // A class that provides the strain rate and its time derivative
    StrainFunction strainFunction;

    // current state variables
    double U, dUdt; // tangential velocity
    double T, dTdt; // temperature
    dvector Y, dYdt; // species mass fractions

    // Diagonalized, linear approximations for terms neglected by splitting
    dvector splitConst; // constant terms
    dvector splitLinear; // diagonal jacobian components

    // Diagonal Jacobian elements for this this term
    dvector diagonalJac;
    bool updateDiagonalJac;

    // Cantera data
    CanteraGas* gas;
    boost::shared_ptr<AdapChem> ckGas;
    bool usingAdapChem;
    perfTimer* reactionRatesTimer;
    perfTimer* thermoTimer;
    perfTimer* jacobianTimer;

    // other parameters
    size_t nSpec;
    int j; // grid index for this system
    dvector W; // species molecular weights [kg/kmol]
    double rhou; // density of the unburned mixture

    // Other quantities
    dvector wDot; // species net production rates [kmol/m^3*s]
    double qDot; // heat release rate per unit volume [W/m^3]

private:
    // Physical properties
    double rho; // density [kg/m^3]
    double cp; // specific heat capacity (average) [J/kg*K]
    dvector cpSpec; // species specific heat capacity [J/mol*K]
    double Wmx; // mixture molecular weight [kg/mol]
    dvector hk; // species enthalpies [J/kmol]

};
