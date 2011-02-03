#pragma once
#include "libconfig.h++"
#include <string>
#include "mathUtils.h"

// indices of the respective equations / solution components
extern const size_t kMomentum, kEnergy, kSpecies, kWmx;

class configOptions
{
public:
    // Populate the members of this class with the contents of the configuration file
    void readOptionsFile(const std::string& filename);

    std::string inputDir;
    std::string outputDir;
    std::string restartFile;
    bool useRelativeRestartPath;

    bool overrideTu;
    bool overrideReactants;
    bool haveRestartFile;

    bool fixedBurnedVal;
    bool fixedLeftLoc;
    bool twinFlame;
    bool curvedFlame;
    bool centeredDifferences;
    bool steadyOnly;
    bool unburnedLeft;

    int regridStepInterval;
    int outputStepInterval;
    int profileStepInterval;
    int currentStateStepInterval; // number of time steps before updating profNow.h5 and outNow.h5
    int terminateStepInterval;
    int integratorRestartInterval;
    double regridTimeInterval;
    double outputTimeInterval;
    double profileTimeInterval;
    double maxTimestep;
    double globalTimestep; // timestep size for the outer solver
    double diffusionTimestep; // timestep for diffusion solver

    double idaRelTol;
    double idaMomentumAbsTol;
    double idaEnergyAbsTol;
    double idaSpeciesAbsTol;

    std::string gasMechanismFile;
    std::string gasPhaseID;
    bool usingMultiTransport;

    // Adapchem configuration
    bool usingAdapChem;
    std::string chemkinMechanismFile;
    std::string adapchemInputFile;
    std::string adapchemModelsFile;
    std::string adapchemDefaultModelFile;
    std::string adapchemDonemodelsFile;
    std::string adapchemRestartFile;
    bool transportEliminationDiffusion;
    bool transportEliminationConvection;
    double adapchem_atol;

    std::string fuel; // molar composition of the fuel
    std::string oxidizer; // molar composition of the oxidizer
    double equivalenceRatio;
    double pressure;

    double Tu;
    double strainRateInitial, strainRateFinal;
    double strainRateDt, strainRateT0;

    // Definition of the starting grid: nPoints evenly spaced between xLeft and xRight
    int nPoints;
    double xLeft, xRight;

    // Tolerances for adaptation and regridding
    double vtol, dvtol, rmTol, dampConst, gridMin, gridMax;
    double uniformityTol, absvtol;
    double boundaryTol, boundaryTolRm;
    int addPointCount;

    double tStart, tEnd;
    bool haveTStart;

    int gridAlpha; // 1 for curved flames, 0 for planar flames

    // Controls which variables are included in the outXXXXXX.h5 files
    bool outputAuxiliaryVariables;
    bool outputTimeDerivatives;
    bool outputHeatReleaseRate;
    bool outputSplitHeatReleaseRate;
    bool outputResidualComponents;
    bool outputProfiles;

    bool terminateForSteadyQdot; // if true, code finishes when integral heat release rate is constant
    double terminationTolerance; // relative tolerance for termination
    double terminationAbsTol; // absolute tolerance
    double terminationPeriod; // period over which to require constant heat release rate
    double terminationMaxTime; // stop integration at this time regardless of the heat release rate

    int errorStopCount;
    bool stopIfError;

    dvector strainRateList;
    bool multiRun; // true if strainRateList is supplied

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

private:
    // load the config option "name" into "value", with "defaultVal" as the default.
    // Return true if a non-default value was found
    template <class T1, class T2>
    bool readOption(const std::string name, T1& value, const T2 defaultVal);

    // Same as readOption, but does not print if the default was used
    template <class T1, class T2>
    bool readOptionQuietDefault(const std::string name, T1& value, const T2 defaultVal);

    libconfig::Config* theConfig;

};
