#pragma once

#include "mathUtils.h"
#include "sundialsUtils.h"
#include "strainFunction.h"
#include "grid.h"
#include "chemistry0d.h"
#include "perfTimer.h"

class ConvectionSystem  : public sdODE, public GridBased
{
    // This is the system representing convection of all state variables in the domain.
public:
    ConvectionSystem();

    // The ODE function: ydot = f(t,y)
    int f(const realtype t, const sdVector& y, sdVector& ydot);
    // This uses an explicit integrator, so the full Jacobian is not needed. However,
    // we still need the diagonal elements to implement the splitting method
    void get_diagonal(const realtype t, dvector& dU, dvector& dT, Array2D& dY);

    int bandedJacobian(const realtype t, const sdVector& y, const sdVector& ydot, sdBandMatrix& J);

    void unroll_y(const sdVector& y); // fill in current state variables from sdvector
    void roll_y(sdVector& y) const; // fill in sdvector with current state variables
    void roll_ydot(sdVector& ydot) const; // fill in sdvector with current time derivatives

    void resize(const size_t nSpec, const size_t nPoints);
    void initialize();

    dvector U, dUdt;
    dvector T, dTdt;
    Array2D Y, dYdt;

    double Tleft; // Temperature left boundary value
    dvector Yleft; // Mass fraction left boundary values
    double rVzero; // mass flux boundary value at j=0

    // Diagonalized, linear approximations for terms neglected by splitting
    dvector splitConstU;
    dvector splitConstT;
    Array2D splitConstY;
    dvector splitLinearU;
    dvector splitLinearT;
    Array2D splitLinearY;

    // Temporaries for the neglected terms
    dvector Tconst;
    dvector Uconst;
    Array2D Yconst;

    // Cantera data
    CanteraGas* gas;

    // Auxiliary variables
    dvector V; // mass flux [kg/m^2*s]
    dvector rV; // (radial) mass flux (r*V) [kg/m^2*s or kg/m*rad*s]

    perfTimer* thermoTimer;

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

class ConvectionSystemUTW : public sdODE, public GridBased
{
    // System representing the coupled convection equations for U, T, and Wmx
    // (tangential velocity, temperature, and mixture molecular weight)
public:
    ConvectionSystemUTW();

    // The ODE function: ydot = f(t,y)
    int f(const realtype t, const sdVector& y, sdVector& ydot);
    // This uses an explicit integrator, so the full Jacobian is not needed. However,
    // we still need the diagonal elements to implement the splitting method
    void get_diagonal(const realtype t, dvector& dU, dvector& dT, dvector& dW);
    int bandedJacobian(const realtype t, const sdVector& y, const sdVector& ydot, sdBandMatrix& J);

    void unroll_y(const sdVector& y); // fill in current state variables from sdvector
    void roll_y(sdVector& y) const; // fill in sdvector with current state variables
    void roll_ydot(sdVector& ydot) const; // fill in sdvector with current time derivatives

    void resize(const size_t nPoints);
    void initialize();

    dvector U, dUdt;
    dvector T, dTdt;
    dvector Wmx, dWdt;

    double Tleft; // Temperature left boundary value
    double Wleft; // mixture molecular weight left boundary value
    double rVzero; // mass flux boundary value at j=0

    // Diagonalized, linear approximations for terms neglected by splitting
    dvector splitConstU;
    dvector splitConstT;
    dvector splitConstW;
    dvector splitLinearU;
    dvector splitLinearT;

    // Temporaries for the neglected terms
    dvector Tconst;
    dvector Uconst;

    // Cantera data
    CanteraGas* gas;

    // Auxiliary variables
    dvector V; // mass flux [kg/m^2*s]
    dvector rV; // (radial) mass flux (r*V) [kg/m^2*s or kg/m*rad*s]

    perfTimer* thermoTimer;

private:
    void V2rV();
    void rV2V();

    size_t nVars; // == 3

    dvector rho; // mixture density [kg/m^3]
    dvector W; // species molecular weights [kg/kmol]

    // variables used internally
    dvector dUdx;
    dvector dTdx;
    dvector dWdx;
};


class ConvectionSystemY : public sdODE, public GridBased
{
    // System representing the convection equation for a single species
    // with a prescribed velocity field V
public:
    ConvectionSystemY();

    // The ODE function: ydot = f(t,y)
    int f(const realtype t, const sdVector& y, sdVector& ydot);
    void get_diagonal(const realtype t, dvector& dU, dvector& dT, Array2D& dY);
    int bandedJacobian(const realtype t, const sdVector& y, const sdVector& ydot, sdBandMatrix& J);

    void resize(const size_t nPoints);
    void initialize();

    // Diagonalized, linear approximations for terms neglected by splitting
    dvector splitConstY;
    dvector splitLinearY;

    // Temporaries for the neglected terms
    dvector Yconst;

    // Auxiliary variables
    dvector v; // velocity [m/s]

private:
    dvector dYdx;
};
