#pragma once

#include <iostream>
#include "sundialsUtils.h"
#include "chemistry0d.h"
#include "grid.h"
#include "readConfig.h"
#include "perfTimer.h"
#include "integrator.h"
#include "sourceSystem.h"
#include "diffusionSystem.h"
#include "convectionSystem.h"

using Cantera::Array2D;
using std::string;

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

    int nPoints;
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

