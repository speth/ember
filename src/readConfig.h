#pragma once
#include "mathUtils.h"

#include <string>

// indices of the respective equations / solution components
const size_t kMomentum = 0;
const size_t kEnergy = 1;
const size_t kSpecies = 2;
const size_t kWmx = 2; // never used in the same systems as kSpecies

//! Possible boundary conditions for the continuity equations
namespace ContinuityBoundaryCondition {
    enum BC {
        Left, //!< V is computed from j = 0
        Right, //!< V is computed from j = jj
        Zero, //!< V = 0 at x = 0
        Qdot, //!< V is calculated from the location of maximum heat release rate
        Temp //!< V is calculated from the temperature midpoint
    };
}

//! Configuration options for the simulation.
/*!
 *  These are more fully documented in `input.py`. For variables that
 *  correspond to a definition in `input.py`, the name used for the variable
 *  in that file is given in square brackets.
 */
class ConfigOptions
{
public:
    ConfigOptions() {}

    //! Returns true if integrator stage data should be saved at the current
    //! time step.
    bool debugIntegratorStages(double t) const;

    //! Set the value of #continuityBC by name
    void setContinuityBC(const std::string& condition);

    std::string outputDir; //!< [paths.outputDir]

    //! [general.fixedBurnedVal] True if burned gas state is fixed at
    //! equilibrium conditions.
    bool fixedBurnedVal;

    //! [general.fixedLeftLocation] True if the position of the leftmost grid
    //! point is fixed.
    bool fixedLeftLoc;

    //! [general.twinFlame] Set to true for the twin flame configuration.
    //! If specified, the leftmost grid point will not extend beyond x = 0.
    bool twinFlame;

    ContinuityBoundaryCondition::BC continuityBC;

    //! [general.curvedFlame] Set to true for the twin flame configuration.
    //! If specified, the leftmost grid point will not extend beyond x = 0.
    //! Also sets #gridAlpha = 1. Otherwise, #gridAlpha = 0.
    bool curvedFlame;
    bool unburnedLeft; //!< [general.unburnedLeft]
    bool fuelLeft; //!< [general.fuelLeft]

    std::string flameType; //!< [initialCondition.flameType]

    int regridStepInterval; //!< [times.regridStepInterval]
    int outputStepInterval; //!< [times.outputStepInterval]
    int profileStepInterval; //!< [times.profileStepInterval]
    int currentStateStepInterval; //!< [times.currentStateStepInterval]
    int terminateStepInterval; //!< [times.terminateStepInterval]
    double regridTimeInterval; //!< [times.regridTimeInterval]
    double outputTimeInterval; //!< [times.outputTimeInterval]
    double profileTimeInterval; //!< [times.profileTimeInterval]
    double globalTimestep; //!< [times.globalTimestep]
    double diffusionTimestepMultiplier;  //!< [times.diffusionTimestepMultiplier]
    std::string splittingMethod; //!< [general.splittingMethod]

    std::string chemistryIntegrator; //!< [general.chemistryIntegrator]
    double integratorRelTol; //!< [cvodeTolerances.relativeTolerance]
    double integratorMomentumAbsTol; //!< [cvodeTolerances.momentumAbsTol]
    double integratorEnergyAbsTol; //!< [cvodeTolerances.energyAbsTol]
    double integratorSpeciesAbsTol; //!< [cvodeTolerances.speciesAbsTol]
    double integratorMinTimestep; //!< [cvodeTolerances.minimumTimestep]

    // QSS integrator parameters
    double qss_epsmax; //!< [qssTolerances.epsmax]
    double qss_epsmin; //!< [qssTolerances.epsmin]
    double qss_dtmin; //!< [qssTolerances.dtmin]
    double qss_dtmax; //!< [qssTolerances.dtmax]
    int qss_iterationCount; //!< [qssTolerances.iterationCount]
    double qss_abstol; //!< [qssTolerances.abstol]
    double qss_minval; //!< [qssTolerances.minval]
    bool qss_stabilityCheck; //!< [qssTolerances.stabilityCheck]

    std::string gasMechanismFile; //!< [chemistry.mechanismFile]
    std::string gasPhaseID; //!< [chemistry.phaseID]
    std::string kineticsModel; //!< [chemistry.kineticsModel]
    std::string transportModel; //!< [chemistry.transportModel]
    double transportThreshold; //!< [chemistry.threshold]
    std::string rateMultiplierFunctionType; //!< [see chemistry.rateMultiplierFunction]

    double pressure; //!< [initialConditoin.pressure]

    std::string strainFunctionType;
    double strainRateInitial; //!< [strainParameters.initial]
    double strainRateFinal; //!< [strainParameters.final]
    double strainRateDt; //!< [strainParameters.dt]
    double strainRateT0; //!< [strainParameters.tStart]

    // Initial profiles specified in configuration file
    dvec x_initial; //!< [initialCondition.x]
    dvec T_initial; //!< [initialCondition.T]
    dvec U_initial; //!< [initialCondition.U]
    dvec V_initial; //!< [initialCondition.V]
    dmatrix Y_initial; //!< [initialCondition.Y]

    //! If false, energy, continuity, and momentum equations are integrated
    //! normally. If true, then time is actually a surrogate for a spatial
    //! variable and the temperature and velocity fields are interpolated from
    //! a supplied data file.
    bool quasi2d;

    //! Data file containing 2D interpolation data [general.interpFile]
    std::string interpFile;

    // Wall Flux properties
    bool wallFlux; //<! True if wall flux is applied when x[0] == 0
    double Tinf; //<! Temperature used for computing wall heat flux [wallFlux.Tinf]
    double Kwall; //<! Thermal conductance of the wall at x = 0 [wallFlux.Kwall]

    // Ignition parameters
    double ignition_tStart; //!< [ignition.tStart]
    double ignition_duration; //!< [ignition.duration]
    double ignition_energy; //!< [ignition.energy]
    double ignition_center; //!< [ignition.center]
    double ignition_stddev; //!< [ignition.stddev]

    // Tolerances for adaptation and regridding
    double vtol; //!< [grid.vtol]
    double dvtol; //!< [grid.dvtol]
    double rmTol; //!< [grid.rmTol]
    double dampConst; //!< [grid.dampConst]
    double gridMin; //!< [grid.gridMin]
    double gridMax; //!< [grid.gridMax]
    double uniformityTol; //!< [grid.uniformityTol]
    double absvtol; //!< [grid.absvtol]
    double boundaryTol; //!< [grid.boundaryTol]
    double boundaryTolRm; //!< [grid.boundaryTolRm]
    double unstrainedDownstreamWidth; //!< [grid.unstrainedDownstreamWidth]
    int addPointCount; //!< [grid.addPointCount]

    double tStart; //!< [times.tStart]
    double tEnd; //!< [terminationCondition.tEnd]
    double tEndMin; //!< [terminationCondition.tMin]
    bool haveTStart;

    int gridAlpha; //!< 1 for curved flames, 0 for planar flames

    // Controls which variables are included in the outXXXXXX.h5 files
    bool outputAuxiliaryVariables; //!< [outputFiles.auxiliaryVariables]
    bool outputTimeDerivatives; //!< [outputFiles.timeDerivatives]
    bool outputHeatReleaseRate; //!< [outputFiles.heatReleaseRate]
    bool outputExtraVariables; //!< [outputFiles.extraVariables]
    bool outputProfiles; //!< [outputFiles.saveProfiles]

    //! Grid point index marked for verbose output (cvodeSteps.py)
    //! [debug.sourcePoint]
    int debugSourcePoint;

    //! Integrator time at which to generate verbose output (then terminate)
    //! [debug.sourceTime]
    double debugSourceTime;

    //! Integrator time at which to begin saving output between integrator
    //! stages [debug.startTime]
    double debugStartTime;

    //! Integrator time at which to stop saving output between integrator
    //! stages [debug.stopTime]
    double debugStopTime;
    bool outputDebugIntegratorStages; //!< [outputFiles.debugIntegratorStages]

    //! if true, code finishes when integral heat release rate is constant
    bool terminateForSteadyQdot;

    //! relative tolerance for termination [terminationCondition.tolerance]
    double terminationTolerance;

    //! absolute tolerance [terminationCondition.abstol]
    double terminationAbsTol;

    //! period over which to require constant heat release rate
    //! [terminationCondition.steadyPeriod]
    double terminationPeriod;

    int errorStopCount; //!< [general.errorStopCount]
    bool stopIfError;
    int nThreads; //!< [general.nThreads]

    //! number of output files written. Initialized with
    //! [outputFiles.firstFileNumber].
    int outputFileNumber;

    //! true if [outputFiles.firstFileNumber] was given in the input file
    bool fileNumberOverride;

    double centerGridMin; //!< [grid.centerGridMin]

    bool xFlameControl;
    double xFlameInitial; //!< [positionControl.xInitial]
    double xFlameFinal; //!< [positionControl.xFinal]
    double xFlameDt; //!< [positionControl.dt]
    double xFlameT0; //!< [positionControl.tStart]
    double xFlameIntegralGain; //!< [positionControl.integralGain]
    double xFlameProportionalGain; //!< [positionControl.proportionalGain]

    bool xStagControl;
    double xStag;
};
