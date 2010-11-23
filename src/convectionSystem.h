#pragma once

#include "mathUtils.h"
#include "sundialsUtils.h"
#include "strainFunction.h"
#include "grid.h"
#include "chemistry0d.h"
#include "perfTimer.h"

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

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
    void get_diagonal(const realtype t, dvector& dU, dvector& dT);
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
    dvector rho; // mixture density [kg/m^3]

    // variables used internally
    dvector dUdx;
    dvector dTdx;
    dvector dWdx;

private:
    void V2rV();
    void rV2V();

    size_t nVars; // == 3
};

typedef std::map<double, dvector> vecInterpolator;

class ConvectionSystemY : public sdODE, public GridBased
{
    // System representing the convection equation for a single species
    // with a prescribed velocity field V
public:
    ConvectionSystemY();

    // The ODE function: ydot = f(t,y)
    int f(const realtype t, const sdVector& y, sdVector& ydot);
    void get_diagonal(const realtype t, dvector& dy);
    int bandedJacobian(const realtype t, const sdVector& y, const sdVector& ydot, sdBandMatrix& J);

    void resize(const size_t nPoints);
    void initialize();

    double Yleft;

    // Diagonalized, linear approximations for terms neglected by splitting
    dvector splitConstY;
    dvector splitLinearY;

    size_t startIndex;
    size_t stopIndex;

    boost::shared_ptr<vecInterpolator> vInterp; // axial (normal) velocity [m/s] at various times

private:
    void update_v(const double t);
    dvector v;
};


class ConvectionSystemSplit : public GridBased
{
    // System which combines a ConvectionSystemUTW objects and several
    // ConvectionSystemY objects that together represent the complete
    // convection term for all components.
public:
    ConvectionSystemSplit();

    void get_diagonal(const realtype t, dvector& dU, dvector& dT, Array2D& dY);

    void setGrid(const oneDimGrid& grid);
    void setTolerances(const configOptions& options);
    void setGas(CanteraGas& gas);
    void resize(const size_t nPointsUTW, const vector<size_t>& nPointsSpec, const size_t nSpec);
    void setSpeciesDomains(vector<size_t>& startIndices, vector<size_t>& stopIndices);
    void setState(const dvector& U, const dvector& T, Array2D& Y);
    void setLeftBC(const double Tleft, const dvector& Yleft);
    void set_rVzero(const double rVzero);
    void initialize(const double t0);
    void evaluate(); // evaluate time derivatives and mass flux at the current state

    // Diagonalized, linear approximations for terms neglected by splitting
    // 'offset' indicates whether constant the values have been offset by the linear term times the current value
    void setSplitConst(const dvector& constU, const dvector& constT, const Array2D& constY, bool offset);
    void setSplitLinear(const dvector& linearU, const dvector& linearT, const Array2D& linearY);
    void setSplitConstU(const dvector& constU);
    void setSplitConstT(const dvector& constT);
    void setSplitConstY(const Array2D& constY, bool offset);
    void setSplitLinearU(const dvector& linearU);
    void setSplitLinearT(const dvector& linearT);
    void setSplitLinearY(const Array2D& linearY);

    void integrateToTime(const double tf);
    void unroll_y(); // convert the solver's solution vectors to the full U, Y, and T
    int getNumSteps();

    dvector U;
    dvector T;
    dvector Wmx;
    Array2D Y;

    // Time derivatives and mass flux are updated by evaluate()
    dvector V;
    dvector dUdt;
    dvector dTdt;
    dvector dWdt;
    Array2D dYdt;

    boost::shared_ptr<vecInterpolator> vInterp;

private:
    // set parameters of a new species solver
    void configureSolver(sundialsCVODE& solver, const size_t k);

    // CVODE integration tolerances
    double reltol; // relative tolerance
    double abstolU; // velocity absolute tolerance
    double abstolT; // temperature absolute tolerance
    double abstolW; // molecular weight absolute tolerance
    double abstolY; // mass fraction absolute tolerance

    ConvectionSystemUTW utwSystem;
    boost::shared_ptr<sundialsCVODE> utwSolver;
    boost::ptr_vector<ConvectionSystemY> speciesSystems;
    boost::ptr_vector<sundialsCVODE> speciesSolvers;

    dvector Yleft;
    dvector W;

    size_t nSpec;
    size_t nVars;
    size_t nPointsUTW;
    vector<size_t> nPointsSpec;

    vector<size_t>* startIndices; // index of leftmost grid point for each component (U,T,Yk)
    vector<size_t>* stopIndices; // index of rightmost grid point for each component (U,T,Yk)

    CanteraGas* gas;
};
