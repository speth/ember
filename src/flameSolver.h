#pragma once

#include <iostream>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>

#include "../adapchem/ckcompat.h"

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

//! Class which manages the main integration loop.
//! Contains the split solvers and is responsible for the large-scale time integration.
class FlameSolver : public GridBased
{
public:
    FlameSolver();
    FlameSolver(const boost::python::api::object& config);
    ~FlameSolver() {}

    void setOptions(const configOptions& options); //!< Set options read from the configuration file
    void initialize(void); //!< call to generate profiles and perform one-time setup
    void run(void); //!< Start the time integration
    void tryrun(void);

    // Called internally by initialize()
    virtual void generateProfile(); //!< Generate initial conditions based on boundary conditions
    virtual void loadProfile(); //!< Load initial conditions from a file

    //! Calculate the mole fraction vector of the reactants based on the
    //! equivalence ratio and the fuel and oxidizer compositions.
    dvector calculateReactantMixture(void);

    //! Determine whether to terminate integration based on reaching steady-state solution.
    bool checkTerminationCondition(void);

    //! Create prof.h5 file.
    //! @param fileName The name of the output file.
    //! @param errorFile Setting this to 'true' will generate a more verbose output file
    //! @param updateDerivatives Setting this to 'true' will compute time derivatives based
    //!        on the current state. This flag should only set if all of the solvers and auxillary
    //!        arrays are sized corresponding to the current state, e.g. after integration but
    //!        but before regridding.
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
    dvector timeVector; //!< Time [s] at the end of every configOptions::outputStepInterval timesteps.
    dvector timestepVector; //!< Global integrator timesteps [s] used to get to the times in #timeVector.
    dvector heatReleaseRate; //!< Integral heat release rate [W/m^2] at the times in #timeVector.
    dvector consumptionSpeed; //!< Consumption speed [m/s] at the times in #timeVector.
    dvector flamePosition; //!< Heat release centroid [m] at the times in #timeVector.

    dvector qDotProd1Vec; //!< Integral heat release rate after the first production half-step.
    dvector qDotProd2Vec; //!< Integral heat release rate after the second production half-step.
    dvector qDotDiff1Vec; //!< Integral heat release rate after the first diffusion half-step.
    dvector qDotDiff2Vec; //!< Integral heat release rate after the second diffusion half-step.
    dvector qDotConv1Vec; //!< Integral heat release rate after the first convection half-step.
    dvector qDotConv2Vec; //!< Integral heat release rate after the second convection half-step.

    configOptions options; //!< Options read from the input file

    double tStart; //!< Integrator start time
    double tEnd; //!< Integrator termination time (upper bound)
    double tNow; //!< Current time reached by the integrator

private:
    boost::ptr_vector<SourceSystem> sourceTerms; // One for each grid point
    boost::ptr_vector<DiffusionSystem> diffusionTerms; // One for each species (plus T and U)

    boost::ptr_vector<sundialsCVODE> sourceSolvers; // One for each grid point
    boost::ptr_vector<BDFIntegrator> diffusionSolvers; // One for each state variable

    boost::ptr_vector<SourceSystemQSS> sourceTermsQSS; // One for each grid point
    vector<bool> useCVODE; // At each grid point, a flag for which integrator to use

    //! System representing the convection terms and containing their integrators
    ConvectionSystemSplit convectionSystem;

    // Boundary values
    double rhou, rhob, rhoLeft, rhoRight;
    double Tu, Tb, Tleft, Tright;
    dvector Yu, Yb, Yleft, Yright;

    void resizeAuxiliary(); //!< Handle resizing of data structures as grid size changes
    void updateCrossTerms(); //!< calculates values of cross-component terms: jSoret, sumcpj, and jCorr
    void updateChemicalProperties(); //!< Update thermodynamic, transport, and kinetic properties
    void updateBC(); //!< Set boundary condition for left edge of domain
    void calculateQdot(); //!< Compute heat release rate using the current temperature and mass fractions
    void calculateTimeDerivatives_old(); //!< Combine dUdt, dTdt and dYdt from the split solvers

    //! Correct the drift of the total mass fractions and reset any negative mass fractions.
    void correctMassFractions();

    //! Calculate the apparent time derivatives for each term in the governing equations
    void calculateTimeDerivatives(double dt);

    // Steps in the Strang split integration process

    // Prepare the balanced splitting terms for integration
    void prepareDiffusionTerms();
    void prepareProductionTerms();
    void prepareConvectionTerms();

    void setDiffusionSolverState(double tInitial);
    void setConvectionSolverState(double tInitial, int stage);
    void setProductionSolverState(double tInitial);

    void extractConvectionState(double dt, int stage);
    void extractDiffusionState(double dt, int stage);
    void extractProductionState(double dt);

    void integrateConvectionTerms(double t, int stage);
    void integrateProductionTerms(double t, int stage);
    void integrateDiffusionTerms(double t, int stage);

    size_t nSpec; //!< Number of chemical species
    size_t nVars; //!< Number of state variables at each grid point (nSpec + 2)
    size_t N; //!< total problem size (nSpec * nVars)

    // State variables:
    dvector U; //!< normalized tangential velocity (u*a/u_inf) [1/s]
    dvector T; //!< temperature [K]
    Array2D Y; //!< species mass fractions, Y(k,j) [-]

    // components of the time derivatives
    Array2D dYdtDiff, dYdtProd, dYdtConv;
    dvector dTdtDiff, dTdtProd, dTdtConv;
    dvector dUdtDiff, dUdtProd, dUdtConv;
    dvector drhodt;

    // State variables at the beginning of the current integrator stage
    dvector Ustart;
    dvector Tstart;
    Array2D Ystart;

    // Changes in each state variable for each terms of the governing equations
    dvector deltaUconv;
    dvector deltaUdiff;
    dvector deltaUprod;
    dvector deltaTconv;
    dvector deltaTdiff;
    dvector deltaTprod;
    Array2D deltaYconv;
    Array2D deltaYdiff;
    Array2D deltaYprod;

    // Auxiliary variables:
    dvector rho; //!< density [kg/m^3]
    dvector jCorr; //!< Correction to ensure sum of mass fractions = 1
    dvector sumcpj; //!< part of the enthalpy flux term
    dvector qDot; //!< Heat release rate [W/m^3]
    Array2D wDot; //!< species production rates [kmol/m^3*s]
    dvector Wmx; //!< mixture molecular weight [kg/kmol]
    dvector W; //!< species molecular weights [kg/kmol]
    dvector mu; //!< dynamic viscosity [Pa*s]
    dvector lambda; //!< thermal conductivity [W/m*K]
    dvector cp; //!< mixture heat capacity [J/kg*K]
    Array2D cpSpec; //!< species molar heat capacities [J/kmol*K]
    Array2D rhoD; //!< density * diffusivity [kg/m*s]
    Array2D Dkt; //!< thermal diffusivity
    Array2D hk; //!< species molar enthalpies [J/kmol]
    Array2D jFick; //!< Fickian mass flux [kg/m^2*s]
    Array2D jSoret; //!< Soret mass flux [kg/m^2*s]

    Array2D dYdtCross; //!< dYdt due to gradients in other species and temperature
    dvector dTdtCross; //!< dTdt due to gradients in gas composition

    // jCorr is a correction to force the net diffusion mass flux to be zero
    // jCorrSystem / jCorrSolver are used to introduce numerical diffusion into
    // jCorr to eliminate spatial instabilities
    DiffusionSystem jCorrSystem;
    BDFIntegrator jCorrSolver;

    // Solver used to determine the subdomain on which to evaluate the
    // transport term for each species
    DiffusionSystem diffusionTestTerm;
    BDFIntegrator diffusionTestSolver;
    vector<size_t> diffusionStartIndices; //!< index of leftmost grid point for each component
    vector<size_t> diffusionStopIndices; //!< index of rightmost grid point for each component
    vector<size_t> nPointsDiffusion; //!< number of grid points for transport of each component

    vector<size_t> convectionStartIndices; //!< index of leftmost grid point for each component
    vector<size_t> convectionStopIndices; //!< index of rightmost grid point for each component
    vector<size_t> nPointsConvection; //!< number of grid points for transport of each component

    //! Function which describes strain rate a(t) and its derivative
    StrainFunction strainfunc;

    double rVcenter; //!< mass flux at centerline [kg/m^2 or kg/m*rad*s]
    double rVzero; //!< mass flux at j=0
    double tFlamePrev, tFlameNext;
    double xFlameTarget, xFlameActual;
    double flamePosIntegralError;

    //! Cantera data
    CanteraGas gas;

    //! Chemkin/AdapChem state mechanism data
    boost::shared_ptr<AdapChem> ckGas;

    void rollVectorVector(vector<dvector>& vv, const dvector& u, const dvector& t, const Array2D& y) const;
    void unrollVectorVector(const vector<dvector>& vv, dvector& u, dvector& t, Array2D& y, size_t i) const;

    void update_xStag(const double t, const bool updateIntError);
    double targetFlamePosition(double t); //!< [m]

    void printPerformanceStats(void);
    void printPerfString(const std::string& label, const perfTimer& T);

    void updateTransportDomain();

    int alpha; //!< curved grid exponent. alpha = 1 for curved flames, 0 for planar flames.
    std::ofstream statsFile;

    // Performance Timers
    // Just the total time:
    perfTimer totalTimer;

    // These add up to the total run time:
    perfTimer setupTimer, splitTimer, reactionTimer, diffusionTimer,
              convectionTimer, regridTimer;

    // These account for special parts of the code
    perfTimer reactionRatesTimer, transportTimer, thermoTimer;
    perfTimer jacobianTimer, adaptiveTransportTimer;
    perfTimer conductivityTimer, viscosityTimer, diffusivityTimer;
};
