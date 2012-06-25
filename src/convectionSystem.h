#pragma once

#include "mathUtils.h"
#include "sundialsUtils.h"
#include "strainFunction.h"
#include "grid.h"
#include "chemistry0d.h"
#include "perfTimer.h"
#include "quasi2d.h"

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

class ConvectionSystemUTW : public sdODE, public GridBased
{
    // System representing the coupled convection equations for U, T, and Wmx
    // (tangential velocity, temperature, and mixture molecular weight)
public:
    ConvectionSystemUTW();

    // The ODE function: ydot = f(t,y)
    int f(const realtype t, const sdVector& y, sdVector& ydot);

    void unroll_y(const sdVector& y); // fill in current state variables from sdvector
    void roll_y(sdVector& y) const; // fill in sdvector with current state variables
    void roll_ydot(sdVector& ydot) const; // fill in sdvector with current time derivatives

    void resize(const size_t nPoints);
    void resetSplitConstants();

    void updateContinuityBoundaryCondition(const dvec& qdot,
                                           ContinuityBoundaryCondition::BC newBC);

    dvec U, dUdt;
    dvec T, dTdt;
    dvec Wmx, dWdt;

    double Tleft; // Temperature left boundary value
    double Wleft; // mixture molecular weight left boundary value
    double rVzero; // mass flux boundary value at j=0

    dvec drhodt;

    // Constant terms introduced by the splitting method
    dvec splitConstT;
    dvec splitConstW;
    dvec splitConstU;

    // Cantera data
    CanteraGas* gas;

    // Auxiliary variables
    dvec V; // mass flux [kg/m^2*s]
    dvec rV; // (radial) mass flux (r*V) [kg/m^2*s or kg/m*rad*s]
    dvec rho; // mixture density [kg/m^3]

    // variables used internally
    dvec dUdx;
    dvec dTdx;
    dvec dWdx;

    ContinuityBoundaryCondition::BC continuityBC;
    size_t jContBC; // The point at which the continuity equation BC is applied
    double xVzero; // Location of the stagnation point (if using fixedZero BC)

private:
    void V2rV();
    void rV2V();

    size_t nVars; // == 3
};

typedef std::map<double, dvec> vecInterpolator;

class ConvectionSystemY : public sdODE, public GridBased
{
    // System representing the convection equation for a single species
    // with a prescribed velocity field V
public:
    ConvectionSystemY() : quasi2d(false) {}

    // The ODE function: ydot = f(t,y)
    int f(const realtype t, const sdVector& y, sdVector& ydot);

    void resize(const size_t nPoints);
    void resetSplitConstants();

    double Yleft;
    int k; // species index (only needed for debugging)
    dvec splitConst; // constant term introduced by splitting

    boost::shared_ptr<vecInterpolator> vInterp; // axial (normal) velocity [m/s] at various times

    //! Interpolators for the quasi-2d problem
    boost::shared_ptr<BilinearInterpolator> vzInterp;
    boost::shared_ptr<BilinearInterpolator> vrInterp;
    bool quasi2d;

private:
    void update_v(const double t);
    dvec v;
};


class ConvectionSystemSplit : public GridBased
{
    // System which combines a ConvectionSystemUTW objects and several
    // ConvectionSystemY objects that together represent the complete
    // convection term for all components.
public:
    ConvectionSystemSplit();

    void setGrid(const oneDimGrid& grid);
    void setTolerances(const configOptions& options);
    void setGas(CanteraGas& gas);
    void resize(const size_t nPoints, const size_t nSpec);
    void setState(const dvec& U, const dvec& T, dmatrix& Y, double tInitial);
    void setLeftBC(const double Tleft, const dvec& Yleft);
    void set_rVzero(const double rVzero);
    void evaluate(); // evaluate time derivatives and mass flux at the current state

    // Time derivatives of species and temperature from the other split terms are needed
    // to correctly compute the density derivative appearing in the continuity equation
    void setDensityDerivative(const dvec& drhodt);

    // Constants introduced by the splitting method
    void setSplitConstants(const dvec& splitConstU,
                           const dvec& splitConstT,
                           const dmatrix& splitConstY);

    void resetSplitConstants();

    void integrateToTime(const double tf);
    void unroll_y(); // convert the solver's solution vectors to the full U, Y, and T
    int getNumSteps();
    void setupQuasi2D(boost::shared_ptr<BilinearInterpolator>& vzInterp,
                      boost::shared_ptr<BilinearInterpolator>& vrInterp);

    dvec U;
    dvec T;
    dvec Wmx;
    dmatrix Y;

    // Time derivatives and mass flux are updated by evaluate()
    dvec V;
    dvec dUdt;
    dvec dTdt;
    dvec dWdt;
    dmatrix dYdt;

    ConvectionSystemUTW utwSystem;

    boost::shared_ptr<vecInterpolator> vInterp;

    PerfTimer utwTimer, speciesTimer;

private:
    // set parameters of a new species solver
    void configureSolver(sundialsCVODE& solver, const size_t k);

    // CVODE integration tolerances
    double reltol; // relative tolerance
    double abstolU; // velocity absolute tolerance
    double abstolT; // temperature absolute tolerance
    double abstolW; // molecular weight absolute tolerance
    double abstolY; // mass fraction absolute tolerance

    boost::shared_ptr<sundialsCVODE> utwSolver;
    boost::ptr_vector<ConvectionSystemY> speciesSystems;
    boost::ptr_vector<sundialsCVODE> speciesSolvers;

    dvec Yleft;
    dvec W;

    size_t nSpec;
    size_t nVars;

    CanteraGas* gas;

    bool quasi2d;
};
