#pragma once

#include <iostream>
#include <boost/ptr_container/ptr_vector.hpp>

#include "readConfig.h"
#include "sundialsUtils.h"
#include "chemistry0d.h"
#include "grid.h"
#include "readConfig.h"
#include "perfTimer.h"
#include "integrator.h"
#include "sourceSystem.h"
#include "diffusionSystem.h"
#include "convectionSystem.h"
#include "strainFunction.h"

using Cantera::Array2D;
using std::string;

class FlameSolver : public GridBased
{
    // This is the system which contains the split solvers and is responsible
    // for the large-scale time integration.
public:
    FlameSolver();
    ~FlameSolver();

    void setOptions(const configOptions& options);
    // Call one of these
    virtual void generateProfile();
    virtual void loadProfile();

    void initialize(void);
    void run(void);

    // Calculate the mole fraction vector of the reactants based on the
    // equivalence ratio and the fuel and oxidizer compositions
    dvector calculateReactantMixture(void);
    bool checkTerminationCondition(void);

    void writeStateFile(const std::string fileName="", bool errorFile=false);

    // Functions for calculating flame information
    double getHeatReleaseRate(void); // [W/m^2]
    double getConsumptionSpeed(void); // [m/s]
    double getFlamePosition(void); // [m]

    // Time-series data
    dvector timeVector; // [s]
    dvector timestepVector; // [s]
    dvector heatReleaseRate; // [W/m^2]
    dvector consumptionSpeed; // [m/s]
    dvector flamePosition; // [m]

    // Options read from the input file
    configOptions options;

    double tStart;
    double tEnd;
    double tNow;

private:
    boost::ptr_vector<SourceSystem> sourceTerms; // One for each grid point
    boost::ptr_vector<DiffusionSystem> diffusionTerms; // One for each species (plus T and U)
    ConvectionSystem convectionTerm;

    boost::ptr_vector<sundialsCVODE> sourceSolvers; // One for each grid point
    boost::ptr_vector<BDFIntegrator> diffusionSolvers; // One for each state variable
    sundialsCVODE* convectionSolver;

    double tPrev;
    double aPrev;

    // Boundary values
    double rhou, rhob, rhoLeft, rhoRight;
    double Tu, Tb, Tleft, Tright;
    dvector Yu, Yb, Yleft, Yright;

    void resizeAuxiliary(); // Handle resizing of data structures as grid size changes
    void updateMinorTerms(); // calculates values of jFick, jSoret, sumcpj, and jCorr
    void updateLeftBC();

    // Utility functions for adaptation & regridding
    void rollVectorVector(const sdVector& y, const dvector& qdot, vector<dvector>& v);
    void unrollVectorVector(const vector<dvector>& v);
    void unrollVectorVectorDot(const vector<dvector>& v);

    // For debugging purposes
    void debugFailedTimestep(const sdVector& y);

    size_t N; // total problem size;
    size_t nSpec; // Number of chemical species
    size_t nVars; // Number of state variables at each grid point (nSpec + 2)

    // State variables:
    dvector U; // normalized tangential velocity (u*a/u_inf) [1/s]
    dvector T; // temperature [K]
    Array2D Y; // species mass fractions, Y(k,j) [-]

    // Time derivatives of state variables:
    dvector dUdt;
    dvector dTdt;
    Array2D dYdt;

    // Extrapolated state from previous timestep:
    dvector Uextrap;
    dvector Textrap;
    Array2D Yextrap;

    // Auxiliary variables:
    dvector rho; // density [kg/m^3]
    dvector jCorr; // Correction to ensure sum of mass fractions = 1
    dvector sumcpj; // part of the enthalpy flux term
    dvector qDot; // Heat release rate [W/m^3]
    Array2D wDot; // species production rates [mol/m^3*s]
    dvector Wmx; // mixture molecular weight
    dvector W;
    dvector mu;
    dvector lambda;
    dvector cp;
    Array2D cpSpec;
    Array2D rhoD;
    Array2D Dkt;
    Array2D hk;
    Array2D jFick;
    Array2D jSoret;
    Array2D jWmx; // diffusion flux proportional to gradient in Wmx

    Array2D constYminor;
    dvector constTminor;

    // Function which describes strain rate a(t) and its derivative
    StrainFunction strainfunc;

    double rVcenter; // mass flux at centerline [kg/m^2 or kg/m*rad*s]
    double rVzero; // mass flux at j=0
    double tFlamePrev, tFlameNext;
    double xFlameTarget, xFlameActual;
    double flamePosIntegralError;

    // Cantera data
    CanteraGas gas;

    void update_xStag(const double t, const bool updateIntError);
    double targetFlamePosition(double t); // [m]

    void V2rV(void);
    void rV2V(void);

    void printPerformanceStats(void);
    void printPerfString(const std::string& label, const perfTimer& T) const;

    // Subdivided governing equation constant terms
    dvector constTdiff, constTconv, constTprod;
    dvector constUdiff, constUconv, constUprod;
    Array2D constYdiff, constYconv, constYprod;

    // Subdivided governing equation linear terms
    dvector linearTdiff, linearTconv, linearTprod;
    dvector linearUdiff, linearUconv, linearUprod;
    Array2D linearYdiff, linearYconv, linearYprod;

    // instantaneous Time derivatives (for output purposes)
    dvector dTdtdiff, dTdtconv, dTdtprod;
    dvector dUdtdiff, dUdtconv, dUdtprod;
    Array2D dYdtdiff, dYdtconv, dYdtprod;

    int alpha;

    double centerVol, centerArea;

    // Performance Timers
    perfTimer perfTimerResFunc, perfTimerPrecondSetup, perfTimerPrecondSolve, perfTimerTransportProps;
    perfTimer perfTimerRxnRates, perfTimerResize, perfTimerLU;
};
