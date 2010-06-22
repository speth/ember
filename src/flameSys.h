#pragma once

#include <iostream>
#include "sundialsUtils.h"
#include "chemistry0d.h"
#include "grid.h"
#include "readConfig.h"
#include "perfTimer.h"
#include "integrator.h"

using Cantera::Array2D;
using std::string;

class StrainFunction
{
public:
    double a(double t) const;
    double dadt(double t) const;
};

class GridBased
{
public:
    GridBased();
    configOptions options;

    // the grid:
    oneDimGrid grid;

private:
    // local names for things that are part of the grid:
    dvector& x;
    dvector& r;
    dvector& rphalf;
    dvector& hh;
    dvector& dlj;
    dvector& cfm;
    dvector& cf;
    dvector& cfp;
};

class SourceSystem : public sdODE
{
    // This is the system representing the (chemical) source term at a point
public:
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

    // other parameters
    size_t nSpec;

private:
    // Cantera data
    CanteraGas* gas;

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

    // The sum of the terms held constant for each component in this system
    dvector constantTerm;
};


class DiffusionSystem : public sdODE, public GridBased
{
    // This is a system representing diffusion of a single solution component
private:
    // Jacobian data
    sdBandMatrix* bandedJacobian;
    vector<int> pMat;

    // Sundials solver parameters
    sdVector* abstol;
    double reltol;

    // The sum of the terms held constant for this system
    dvector constantTerm;
};

class SpeciesDiffusionSystem : public DiffusionSystem
{
    // This system represents the diffusion of a single species at all grid points
public:
    int f(realtype t, sdVector& y, sdVector& ydot);

    int preconditionerSetup(realtype t, sdVector& yIn, sdVector& ydotIn,
                            sdVector& resIn, realtype c_j);

    int preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn, sdVector& resIn,
                            sdVector& rhs, sdVector& outVec, realtype c_j, realtype delta);

    // Diffusion coefficients
    dvector rhoD; // density-weighted, mixture-averaged diffusion coefficients [kg/m*s] (= rho*Dkm)
    dvector Dkt; // thermal diffusion coefficients [kg/m*s]

    // Diffusion mass fluxes
    dvector jFick; // Normal diffusion (Fick's Law) [kg/m^2*s]
    dvector jSoret; // Soret Effect diffusion [kg/m^2*s]
    dvector dYkdx; // upwinded
};


class TemperatureDiffusionSystem : public DiffusionSystem
{
    // This system represents the diffusion Temperature at all grid points
public:
    int f(realtype t, sdVector& y, sdVector& ydot);

    int preconditionerSetup(realtype t, sdVector& yIn, sdVector& ydotIn,
                            sdVector& resIn, realtype c_j);

    int preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn, sdVector& resIn,
                            sdVector& rhs, sdVector& outVec, realtype c_j, realtype delta);

private:
    dvector sumcpj; // for enthalpy flux term [W/m^2*K]
    dvector lambda; // thermal conductivity [W/m*K]
    dvector cp; // specific heat capacity (average) [J/kg*K]
    Array2D cpSpec; // species specific heat capacity [J/mol*K]
    dvector W; // species molecular weights [kg/kmol]
    dvector qFourier; // heat flux [W/m^2]
    dvector dTdx; // upwinded
    dvector dTdxCen; // centered difference (for enthalpy flux term)
};


class MomentumDiffusionSystem : public DiffusionSystem
{
    // This system represents the diffusion of tangential momentum at all grid points
public:
    int f(realtype t, sdVector& y, sdVector& ydot);

    int preconditionerSetup(realtype t, sdVector& yIn, sdVector& ydotIn,
                            sdVector& resIn, realtype c_j);

    int preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn, sdVector& resIn,
                            sdVector& rhs, sdVector& outVec, realtype c_j, realtype delta);


private:
    dvector mu; // viscosity [Pa*s]
    dvector dUdx; // upwinded
};


class ConvectionSystem : public sdODE, public GridBased
{
    // This is the system representing convection of all state variables in the domain.
public:
    int f(realtype t, sdVector& y, sdVector& ydot);
    // This uses an explicit integrator, so no Jacobian/Preconditioner is necessary

    // Sundials solver parameters
    sdVector* abstol;
    double reltol;

private:
    dvector drhodt;
    dvector Wmx; // mixture molecular weight [kg/kmol]
    dvector rV; // (radial) mass flux (r*V) [kg/m^2*s or kg/m*rad*s]
    dvector W; // species molecular weights [kg/kmol]

    // The sum of the terms held constant for this system
    dvector constantTerm;
};

class FlameSystem : public GridBased
{
    // This is the system which contains the split solvers and is responsible
    // for the large-scale time integration.
public:
    FlameSystem();
    ~FlameSystem();

private:
    vector<SourceSystem*> sourceTerms; // One for each grid point
    vector<DiffusionSystem*> diffusionTerms; // One for each species (plus T and U)
    ConvectionSystem convectionTerm;

    // Problem definition
    std::string reactants;
    std::string diluent;
    double xLeft, xRight;
    int nPoints;
    void generateInitialProfiles(void);
    void loadInitialProfiles(void);

    double tStart;
    double tEnd;
    double tNow;
    double tPrev;
    double aPrev;

    // Boundary values
    double rhou, rhob, rhoLeft, rhoRight;
    double Tu, Tb, Tleft, Tright;
    dvector Yu, Yb, Yleft, Yright;

    void setup(void);
    void copyOptions(void);

    // Utility functions for adaptation & regridding
    void rollVectorVector(const sdVector& y, const dvector& qdot, vector<dvector>& v);
    void unrollVectorVector(const vector<dvector>& v);
    void unrollVectorVectorDot(const vector<dvector>& v);

    void printForMatlab(ofstream& file, dvector& v, int index, char* name);
    void writeStateFile(const std::string fileName="", bool errorFile=false);

    // For debugging purposes
    void debugFailedTimestep(const sdVector& y);

    // these should be read-only:
    int N; // total problem size;
    int nVars; // Number of solution variables at each point
    int nSpec; // Number of chemical species

    // State variables:
    dvector U; // normalized tangential velocity (u*a/u_inf) [1/s]
    dvector T; // temperature [K]
    Array2D Y; // species mass fractions, Y(k,j)

    // Time derivatives of state variables:
    dvector dVdt;
    dvector dUdt;
    dvector dTdt;
    Array2D dYdt;

    // Auxiliary variables:
    dvector rho; // density [kg/m^3]
    dvector jCorr; // Correction to ensure sum of mass fractions = 1

    // Strain rate parameters:
    // The strain rate is constant at a=strainRateInitial until
    // t=strainRateT0, after which it increases linearly to a=strainRateFinal
    // at t=strainRateT0+strainRateDt, after which it remains constant.
    double strainRateInitial, strainRateFinal; // [1/s]
    double strainRateDt; // [s]
    double strainRateT0; // [s]

    double rVcenter; // mass flux at centerline [kg/m^2 or kg/m*rad*s]
    double rVzero; // mass flux at j=0
    double tFlamePrev, tFlameNext;
    double xFlameTarget, xFlameActual;
    double flamePosIntegralError;

    // Cantera data
    CanteraGas gas;

    // Miscellaneous options
    configOptions options;

    // Functions for calculating flame information
    double getHeatReleaseRate(void); // [W/m^2]
    double getConsumptionSpeed(void); // [m/s]
    double getFlamePosition(void); // [m]

    double strainRate(const double t); // [1/s]
    double dStrainRatedt(const double t); // [1/s^2]

    void update_xStag(const double t, const bool updateIntError);
    double targetFlamePosition(double t); // [m]

    void V2rV(void);
    void rV2V(void);

    void printPerformanceStats(void);
    void printPerfString(const std::string& label, const perfTimer& T) const;

    // Subdivided governing equation components
    dvector energyDiff, energyConv, energyProd;
    dvector momentumDiff, momentumConv, momentumProd;
    Array2D speciesDiff, speciesConv, speciesProd;

    int alpha;

    double centerVol, centerArea;

    // Performance Timers
    perfTimer perfTimerResFunc, perfTimerPrecondSetup, perfTimerPrecondSolve, perfTimerTransportProps;
    perfTimer perfTimerRxnRates, perfTimerSetup, perfTimerLU;
};


class flameSys
{
public:
    // Remnants of the old class flameSys that should be moved somewhere else or deleted

    // Utility functions
    void unrollY(const sdVector& y);
    void unrollYdot(const sdVector& yDot);
    void rollY(sdVector& y);
    void rollYdot(sdVector& yDot);

    void updateTransportProperties(void);
    void updateThermoProperties(void);
    void updateLeftBC(void);

    dvector V; // mass flux normal to flame per unit area (rho*v) [kg/m^2*s]
    bool inTestPreconditioner;
    void testPreconditioner(void);
};

