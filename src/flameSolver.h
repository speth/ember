#pragma once

#include "readConfig.h"
#include "chemistry0d.h"
#include "grid.h"
#include "splitSolver.h"
#include "readConfig.h"
#include "perfTimer.h"
#include "integrator.h"
#include "sourceSystem.h"
#include "diffusionSystem.h"
#include "convectionSystem.h"
#include "quasi2d.h"

#include <iostream>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/enumerable_thread_specific.h"
#include "tbb/mutex.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/tbb_exception.h"

using std::string;
class SourceTermWrapper;
class DiffusionTermWrapper;
class StrainFunction;

//! Class which manages the main integration loop.
//! Contains the split solvers and is responsible for the large-scale time integration.
class FlameSolver : public GridBased, public SplitSolver
{
public:
    FlameSolver();
    virtual ~FlameSolver();

    void setOptions(const ConfigOptions& options); //!< Set options read from the configuration file
    void initialize(); //!< call to generate profiles and perform one-time setup
    void finalize();

    //! Load initial temperature, mass fraction and velocity profiles.
    //! Profiles are loaded either from an HDF5 "restart" file or from the
    //! the ConfigOptions object, populated from the Python "Config" object.
    virtual void loadProfile();

    //! Determine whether to terminate integration based on reaching steady-state solution.
    bool checkTerminationCondition(void);

    void writeStateFile(const std::string& fileName="", bool errorFile=false, bool updateDerivatives=true);
    void writeTimeseriesFile(const std::string& filename); //!< create out.h5 file

    //! Compute integral heat release rate [W/m^2].
    //! Assumes that #qDot has already been updated.
    double getHeatReleaseRate(void);

    //! Compute flame consumption speed [m/s]. Computed by integrating the energy equation:
    //! \f[ S_c = \frac{\int_{-\infty}^\infty \dot{q}/c_p}{\rho_u(T_b-T_u)} \f]
    //! Assumes that #qDot has already been updated.
    double getConsumptionSpeed(void);

    //! Compute position of flame (heat release centroid) [m]
    double getFlamePosition(void);

    // Time-series data
    dvector timeVector; //!< Time [s] at the end of every ConfigOptions::outputStepInterval timesteps.
    dvector timestepVector; //!< Global integrator timesteps [s] used to get to the times in #timeVector.
    dvector heatReleaseRate; //!< Integral heat release rate [W/m^2] at the times in #timeVector.
    dvector consumptionSpeed; //!< Consumption speed [m/s] at the times in #timeVector.
    dvector flamePosition; //!< Heat release centroid [m] at the times in #timeVector.
    dvector strainRateVector; //!< Strain rate [1/s] at the times in #timeVector.
    dvector dStrainRateVector; //!< Time derivative of strain rate [1/s^2] at the times in #timeVector.

    long int nTotal; //!< total number of timesteps taken
    int nRegrid; //!< number of time steps since regridding/adaptation
    int nOutput; //!< number of time steps since storing integral flame parameters
    int nProfile; //!< number of time steps since saving flame profiles
    int nTerminate; //!< number of steps since last checking termination condition
    int nCurrentState; //!< number of time steps since profNow.h5 and out.h5 were written
    double terminationCondition; //!< Current value to be compared with terminationTolerance

    double tOutput; //!< time of next integral flame parameters output (this step)
    double tRegrid; //!< time of next regridding
    double tProfile; //!< time of next profile output

    boost::ptr_vector<SourceSystem> sourceTerms; // One for each grid point
    vector<bool> useCVODE; // At each grid point, a flag for which integrator to use

    boost::ptr_vector<DiffusionSystem> diffusionTerms; // One for each species (plus T and U)
    boost::ptr_vector<TridiagonalIntegrator> diffusionSolvers; // One for each state variable

    //! System representing the convection terms and containing their integrators
    ConvectionSystemSplit convectionSystem;

    // Boundary values
    double rhou, rhob, rhoLeft, rhoRight;
    double Tu, Tb, Tleft, Tright;
    dvec Yu, Yb, Yleft, Yright;

    void setupStep();
    void prepareIntegrators();
    int finishStep();
    void resizeAuxiliary(); //!< Handle resizing of data structures as grid size changes
    void resizeMappedArrays(); //!< update data that shadows SplitSolver arrays
    void updateCrossTerms(); //!< calculates values of cross-component terms: jSoret, sumcpj, and jCorr
    void updateChemicalProperties(); //!< Update thermodynamic, transport, and kinetic properties
    void updateBC(); //!< Set boundary condition for left edge of domain
    void calculateQdot(); //!< Compute heat release rate using the current temperature and mass fractions

    //! Correct the drift of the total mass fractions and reset any negative mass fractions.
    void correctMassFractions();

    void setDiffusionSolverState(double tInitial);
    void setConvectionSolverState(double tInitial);
    void setProductionSolverState(double tInitial);

    // Steps in the Strang split integration process
    void integrateConvectionTerms(double tStart, double tEnd);
    void integrateProductionTerms(double tStart, double tEnd);
    void integrateDiffusionTerms(double tStart, double tEnd         );

    size_t nSpec; //!< Number of chemical species
    size_t nVars; //!< Number of state variables at each grid point (nSpec + 2)

    // State variables:
    VecMap U; //!< normalized tangential velocity (u*a/u_inf) [1/s]
    VecMap T; //!< temperature [K]
    MatrixMap Y; //!< species mass fractions, Y(k,j) [-]

    // Auxiliary variables:
    dvec rho; //!< density [kg/m^3]
    dvec drhodt; //!< time derivative of density [kg/m^3*s]
    dvec jCorr; //!< Correction to ensure sum of mass fractions = 1
    dvec sumcpj; //!< part of the enthalpy flux term
    dvec qDot; //!< Heat release rate [W/m^3]
    dmatrix wDot; //!< species production rates [kmol/m^3*s]
    dvec Wmx; //!< mixture molecular weight [kg/kmol]
    dvec W; //!< species molecular weights [kg/kmol]
    dvec mu; //!< dynamic viscosity [Pa*s]
    dvec lambda; //!< thermal conductivity [W/m*K]
    dvec cp; //!< mixture heat capacity [J/kg*K]
    dmatrix cpSpec; //!< species molar heat capacities [J/kmol*K]
    dmatrix rhoD; //!< density * diffusivity [kg/m*s]
    dmatrix Dkt; //!< thermal diffusivity
    dmatrix hk; //!< species molar enthalpies [J/kmol]
    dmatrix jFick; //!< Fickian mass flux [kg/m^2*s]
    dmatrix jSoret; //!< Soret mass flux [kg/m^2*s]

    // jCorr is a correction to force the net diffusion mass flux to be zero
    // jCorrSystem / jCorrSolver are used to introduce numerical diffusion into
    // jCorr to eliminate spatial instabilities
    DiffusionSystem jCorrSystem;
    TridiagonalIntegrator jCorrSolver;

    //! Function which describes strain rate a(t) and its derivative
    StrainFunction* strainfunc;

    double rVcenter; //!< mass flux at centerline [kg/m^2 or kg/m*rad*s]
    double rVzero; //!< mass flux at j=0
    double tFlamePrev, tFlameNext;
    double xFlameTarget, xFlameActual;
    double flamePosIntegralError;

    //! Cantera data
    CanteraGas gas;

    tbb::enumerable_thread_specific<CanteraGas> gases;
    tbb::mutex gasInitMutex;
    tbb::task_scheduler_init tbbTaskSched;

    void rollVectorVector(vector<dvector>& vv, const dmatrix& M) const;
    void unrollVectorVector(vector<dvector>& vv, dmatrix& M, size_t i) const;

    void update_xStag(const double t, const bool updateIntError);
    double targetFlamePosition(double t); //!< [m]

    void printPerformanceStats(void);
    void printPerfString(std::ostream& stats, const std::string& label, const PerfTimer& T);

    //! Data for solving quasi-2D method-of-lines problems
    boost::shared_ptr<BilinearInterpolator> vzInterp, vrInterp, TInterp;

    // Performance Timers
    // Just the total time:
    PerfTimer totalTimer;

    // These add up to the total run time:
    PerfTimer setupTimer, reactionTimer, diffusionTimer,
              convectionTimer, regridTimer;

    // These account for special parts of the code
    PerfTimer reactionRatesTimer, transportTimer, thermoTimer;
    PerfTimer jacobianTimer;
    PerfTimer conductivityTimer, viscosityTimer, diffusivityTimer;

    friend class SourceTermWrapper;
    friend class DiffusionTermWrapper;
};

class SourceTermWrapper {
public:
    SourceTermWrapper(FlameSolver* parent_, double t_)
        : parent(parent_)
        , t(t_)
        {}
    void operator()(const tbb::blocked_range<size_t>& r) const;
private:
    FlameSolver* parent;
    int stage;
    double t;
};

class DiffusionTermWrapper {
public:
    DiffusionTermWrapper(FlameSolver* parent_, double t_)
        : parent(parent_)
        , t(t_)
        {}
    void operator()(const tbb::blocked_range<size_t>& r) const;
private:
    FlameSolver* parent;
    double t;
};
