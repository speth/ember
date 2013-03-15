#pragma once
#include "mathUtils.h"

#include <string>

// indices of the respective equations / solution components
extern const size_t kMomentum, kEnergy, kSpecies, kWmx;

namespace ContinuityBoundaryCondition {
    enum BC {
        Left, // V is computed from j = 0
        Right, // V is computed from j = jj
        Zero, // V = 0 at x = 0
        Qdot, // V is calculated from the location of maximum heat release rate
        Temp // V is calculated from the temperature midpoint
    };
}

class ConfigOptions
{
public:
    ConfigOptions() {}

    //! Returns true if integrator stage data should be saved at the current
    //! time step.
    bool debugIntegratorStages(double t) const;

    //! Set the value of #continuityBC by name
    void setContinuityBC(const std::string& condition);

    //! [paths.outputDir] Directory to store output files.
    //! Automatically created if it does not already exist.
    std::string outputDir;

    //! [initialCondition.file] A prof.h5 file to load initial conditions from.
    std::string restartFile;

    //! Set a new temperature boundary value. Automatically set to 'true' if initialCondition.Tu
    //! is specified in the config file.
    bool overrideTu;

    //! Set new species boundary values. Automatically set to 'true' if initialCondition.fuel
    //! is specified in the config file.
    bool overrideReactants;

    bool haveRestartFile; //!< True if a restart file has been specified.

    //! [general.fixedBurnedVal] True if burned gas state is fixed at equilibrium conditions.
    bool fixedBurnedVal;

    //! [general.fixedLeftLocation] True if the position of the leftmost grid point is fixed.
    bool fixedLeftLoc;

    //! [general.twinFlame] Set to true for the twin flame configuration.
    //! If specified, the leftmost grid point will not extend beyond x = 0.
    bool twinFlame;

    ContinuityBoundaryCondition::BC continuityBC;

    //! [general.curvedFlame] Set to true for the twin flame configuration.
    //! If specified, the leftmost grid point will not extend beyond x = 0.
    //! Also sets #gridAlpha = 1. Otherwise, #gridAlpha = 0.
    bool curvedFlame;
    bool unburnedLeft;
    bool fuelLeft;

    std::string flameType; //!< "premixed" or "diffusion"

    int regridStepInterval;
    int outputStepInterval;
    int profileStepInterval;
    int currentStateStepInterval; // number of time steps before updating profNow.h5 and out.h5
    int terminateStepInterval;
    double regridTimeInterval;
    double outputTimeInterval;
    double profileTimeInterval;
    double globalTimestep; // timestep size for the outer solver
    double diffusionTimestepMultiplier; // timestep for diffusion solver
    std::string splittingMethod;

    std::string chemistryIntegrator;
    double integratorRelTol;
    double integratorMomentumAbsTol;
    double integratorEnergyAbsTol;
    double integratorSpeciesAbsTol;
    double integratorMinTimestep; //!< minimum timestep for the chemistry integrator [s]

    // QSS integrator parameters
    double qss_epsmax;
    double qss_epsmin;
    double qss_dtmin;
    double qss_dtmax;
    int qss_iterationCount;
    double qss_abstol;
    double qss_minval;
    bool qss_stabilityCheck;

    std::string gasMechanismFile;
    std::string gasPhaseID;
    std::string transportModel;
    double transportThreshold;

    std::string fuel; // molar composition of the fuel
    std::string oxidizer; // molar composition of the oxidizer
    double equivalenceRatio;
    double pressure;

    double Tu; //!< Temperature of the unburned mixture for premixed flames
    double Tfuel; //!< Temperature of the fuel for diffusion flames
    double Toxidizer; //!< Temperature of the oxidizer for diffusion flames

    std::string strainFunctionType;
    double strainRateInitial, strainRateFinal;
    double strainRateDt, strainRateT0;

    // Initial profiles specified in configuration file
    bool haveInitialProfiles;
    dvec x_initial;
    dvec T_initial;
    dvec U_initial;
    dmatrix Y_initial;
    double rVzero_initial;

    //! If false, energy, continuity, and momentum equations are integrated
    //! normally. If true, then time is actually a surrogate for a spatial
    //! variable and temperature velocity fields are interpolated from a
    //! supplied data file.
    bool quasi2d;
    std::string interpFile; //!< data file containing 2D interpolation data

    // Wall Flux properties
    bool wallFlux; //<! True if wall flux is applied when x[0] == 0
    double Tinf; //<! Temperature used for computing wall heat flux
    double Kwall; //<! Thermal conductance of the wall at x = 0

    // Ignition parameters
    double ignition_tStart;
    double ignition_duration;
    double ignition_energy;
    double ignition_center;
    double ignition_stddev;

    // Tolerances for adaptation and regridding
    double vtol, dvtol, rmTol, dampConst, gridMin, gridMax;
    double uniformityTol, absvtol;
    double boundaryTol, boundaryTolRm;
    double unstrainedDownstreamWidth;
    int addPointCount;

    double tStart, tEnd;
    bool haveTStart;

    int gridAlpha; // 1 for curved flames, 0 for planar flames

    // Controls which variables are included in the outXXXXXX.h5 files
    bool outputAuxiliaryVariables;
    bool outputTimeDerivatives;
    bool outputHeatReleaseRate;
    bool outputExtraVariables;
    bool outputProfiles;

    int debugSourcePoint; //!< Grid point index marked for verbose output (cvodeSteps.py)
    double debugSourceTime; //!< Integrator time at which to generate verbose output (then terminate)
    double debugStartTime; //!< Integrator time at which to begin saving output between integrator stages
    double debugStopTime; //!< Integrator time at which to stop saving output between integrator stages
    bool outputDebugIntegratorStages;

    bool terminateForSteadyQdot; // if true, code finishes when integral heat release rate is constant
    double terminationTolerance; // relative tolerance for termination
    double terminationAbsTol; // absolute tolerance
    double terminationPeriod; // period over which to require constant heat release rate

    int errorStopCount;
    bool stopIfError;
    int nThreads;

    int outputFileNumber; // number of output files written
    bool fileNumberOverride; // true if outputFileNumbe was given in the input file

    double centerGridMin;

    bool xFlameControl;
    double xFlameInitial, xFlameFinal;
    double xFlameDt, xFlameT0;
    double xFlameIntegralGain;
    double xFlameProportionalGain;

    bool xStagControl;
    double xStag;

};
