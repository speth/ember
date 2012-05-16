#include "readConfig.h"
#include "boost/filesystem.hpp"
#include "debugUtils.h"

#include <boost/python.hpp>

const size_t kMomentum = 0;
const size_t kEnergy = 1;
const size_t kSpecies = 2;
const size_t kWmx = 2; // never used in the same systems as kSpecies

namespace bp = boost::python;

template <class T1>
bool configOptions::readOption(const bp::object& conf, const char* name, T1& value)
{
    const bp::object& obj = conf.attr(name);
    if (obj.ptr() != Py_None) {
        value = bp::extract<T1>(obj);
        return true;
    } else {
        return false;
    }
}

// Read 1D array (list, numpy array, etc.) into a vector
template <>
bool configOptions::readOption(const bp::object& conf,
                               const char* name, dvector& value)
{
    const bp::object& obj = conf.attr(name);
    if (obj.ptr() != Py_None) {
        int N = bp::len(obj);
        value.resize(N);
        for (int i=0; i<N; i++) {
            value[i] = bp::extract<double>(obj[i]);
        }
        return true;
    } else {
        return false;
    }
}

// Read 1D array (list, numpy array, etc.) into an  Eigen vector
template <>
bool configOptions::readOption(const bp::object& conf,
                               const char* name, dvec& value)
{
    const bp::object& obj = conf.attr(name);
    if (obj.ptr() != Py_None) {
        int N = bp::len(obj);
        value.resize(N);
        for (int i=0; i<N; i++) {
            value[i] = bp::extract<double>(obj[i]);
        }
        return true;
    } else {
        return false;
    }
}

// Read 2D Numpy array into an Eigen matrix
template <>
bool configOptions::readOption(const bp::object& conf,
                               const char* name, dmatrix& value)
{
    const bp::object& obj = conf.attr(name);
    if (obj.ptr() != Py_None) {
        int N = bp::extract<int>(obj.attr("shape")[0]);
        int M = bp::extract<int>(obj.attr("shape")[1]);
        value.resize(N,M);
        for (int i=0; i<N; i++) {
            for (int j=0; j<M; j++) {
                value(i,j) = bp::extract<double>(obj[bp::make_tuple(i,j)]);
            }
        }
        return true;
    } else {
        return false;
    }
}

configOptions::configOptions(const bp::object& conf)
{
    using bp::object;
    using std::string;

    object None;

    // Paths
    const object& paths = conf.attr("paths");
    readOption(paths, "inputDir", inputDir);
    readOption(paths, "outputDir", outputDir);

    string logFileName;
    if (readOption(paths, "logFile", logFileName)) {
        logFile.open(logFileName);
    }

    // General
    const object& general = conf.attr("general");
    readOption(general, "fixedBurnedVal", fixedBurnedVal);
    readOption(general, "fixedLeftLocation", fixedLeftLoc);
    readOption(general, "unburnedLeft", unburnedLeft);
    readOption(general, "fuelLeft", fuelLeft);
    readOption(general, "nThreads", nThreads);

    readOption(general, "curvedFlame", curvedFlame);
    readOption(general, "twinFlame", twinFlame);

    readOption(general, "interpFile", interpFile);

    stopIfError = readOption(general, "errorStopCount", errorStopCount);
    readOption(general, "chemistryIntegrator", chemistryIntegrator);
    readOption(general, "splittingMethod", splittingMethod);

    std::string bcString;
    readOption(general, "continuityBC", bcString);
    if (bcString == "fixedLeft") {
        continuityBC = ContinuityBoundaryCondition::Left;
    } else if (bcString == "fixedRight") {
        continuityBC = ContinuityBoundaryCondition::Right;
    } else if (bcString == "fixedQdot") {
        continuityBC = ContinuityBoundaryCondition::Qdot;
    } else if (bcString == "fixedTemperature") {
        continuityBC = ContinuityBoundaryCondition::Temp;
    } else if (bcString == "stagnationPoint") {
        continuityBC = ContinuityBoundaryCondition::Zero;
    }

    // Chemistry
    const object& chem = conf.attr("chemistry");
    readOption(chem, "mechanismFile", gasMechanismFile);
    readOption(chem, "phaseID", gasPhaseID);
    readOption(chem, "transportModel", transportModel);
    readOption(chem, "threshold", transportThreshold);

    // Transport species elimination
    const object& transport = conf.attr("transportElimination");
    readOption(transport, "atol", transport_atol);
    readOption(transport, "diffusion", transportEliminationDiffusion);
    readOption(transport, "convection", transportEliminationConvection);
    readOption(transport, "stepInterval", transportEliminationStepInterval);

    // Initial Condition
    const object& ic = conf.attr("initialCondition");
    readOption(ic, "nPoints", nPoints);
    readOption(ic, "xLeft", xLeft);
    readOption(ic, "xRight", xRight);
    readOption(ic, "centerWidth", initialCenterWidth);
    readOption(ic, "slopeWidth", initialSlopeWidth);
    readOption(ic, "smoothCount", initialSmoothCount);
    haveRestartFile = readOption(ic, "restartFile", restartFile);
    useRelativeRestartPath = true;
    readOption(ic, "relativeRestartPath", useRelativeRestartPath);
    if (haveRestartFile) {
        haveRestartFile = boost::filesystem::exists(inputDir + "/" + restartFile);
        if (!haveRestartFile) {
            logFile.write(format("WARNING: couldn't find restart file '%s/%s'") %
                inputDir % restartFile);
        }
    }

    // Initial Condition: full profiles
    bool test = readOption(ic, "haveProfiles", haveInitialProfiles);
    haveInitialProfiles = haveInitialProfiles && test;
    if (haveInitialProfiles) {
        readOption(ic, "x", x_initial);
        readOption(ic, "T", T_initial);
        readOption(ic, "U", U_initial);
        readOption(ic, "Y", Y_initial);
        readOption(ic, "rVzero", rVzero_initial);
    }

    // Wall flux boundary condition
    const object& wall = conf.attr("wallFlux");
    if (wall != None) {
        wallFlux = true;
        readOption(wall, "Tinf", Tinf);
        readOption(wall, "Kwall", Kwall);
    } else {
        wallFlux = false;
        Tinf = Tu;
        Kwall = 0;
    }

    overrideTu = readOption(ic, "Tu", Tu);
    overrideReactants = readOption(ic, "fuel", fuel);
    readOption(ic, "flameType", flameType);
    readOption(ic, "oxidizer", oxidizer);
    readOption(ic, "Tfuel", Tfuel);
    readOption(ic, "Toxidizer", Toxidizer);
    readOption(ic, "equivalenceRatio", equivalenceRatio);
    readOption(ic, "pressure", pressure);
    quasi2d = (flameType == "quasi2d");

    // Ignition parameters
    const object& ignition = conf.attr("ignition");
    readOption(ignition, "tStart", ignition_tStart);
    readOption(ignition, "duration", ignition_duration);
    readOption(ignition, "energy", ignition_energy);
    readOption(ignition, "center", ignition_center);
    readOption(ignition, "stddev", ignition_stddev);

    // Strain Rate Parameters
    const object& strain = conf.attr("strainParameters");
    readOption(strain, "initial", strainRateInitial);
    readOption(strain, "final", strainRateFinal);
    readOption(strain, "tStart", strainRateT0);
    readOption(strain, "dt", strainRateDt);

    const object& positionControl = conf.attr("positionControl");
    if (positionControl != None) {
        xStagControl = true;
    } else {
        xStagControl = false;
    }

    if (positionControl != None &&
        readOption(positionControl, "xInitial", xFlameInitial) &&
        readOption(positionControl, "xFinal", xFlameFinal) &&
        readOption(positionControl, "tStart", xFlameT0) &&
        readOption(positionControl, "dt", xFlameDt)) {
        readOption(positionControl, "integralGain", xFlameIntegralGain);
        readOption(positionControl, "proportionalGain", xFlameProportionalGain);
        xFlameControl = true;
        xStagControl = true;
    } else {
        xFlameControl = false;
    }

    // Grid
    const object& grid = conf.attr("grid");
    readOption(grid, "centerGridMin", centerGridMin);
    readOption(grid, "vtol", vtol);
    readOption(grid, "dvtol", dvtol);
    readOption(grid, "rmTol", rmTol);
    readOption(grid, "dampConst", dampConst);
    readOption(grid, "gridMin", gridMin);
    readOption(grid, "gridMax", gridMax);
    readOption(grid, "uniformityTol", uniformityTol);
    readOption(grid, "absvtol", absvtol);
    readOption(grid, "boundaryTol", boundaryTol);
    readOption(grid, "boundaryTolRm", boundaryTolRm);
    readOption(grid, "addPointCount", addPointCount);
    readOption(grid, "unstrainedDownstreamWidth", unstrainedDownstreamWidth);

    // Times
    const object& times = conf.attr("times");
    haveTStart = readOption(times, "tStart", tStart);
    readOption(times, "regridTimeInterval", regridTimeInterval);
    readOption(times, "regridStepInterval", regridStepInterval);
    readOption(times, "outputTimeInterval", outputTimeInterval);
    readOption(times, "outputStepInterval", outputStepInterval);
    readOption(times, "profileTimeInterval", profileTimeInterval);
    readOption(times, "profileStepInterval", profileStepInterval);
    readOption(times, "currentStateStepInterval", currentStateStepInterval);
    readOption(times, "terminateStepInterval", terminateStepInterval);
    readOption(times, "globalTimestep", globalTimestep);
    readOption(times, "diffusionTimestepMultiplier", diffusionTimestepMultiplier);

    // Debug parameters
    const object& debug = conf.attr("debug");
    readOption(debug, "adaptation", debugParameters::debugAdapt);
    readOption(debug, "regridding", debugParameters::debugRegrid);
    readOption(debug, "timesteps", debugParameters::debugTimesteps);
    readOption(debug, "flameRadiusControl", debugParameters::debugFlameRadiusControl);
    #ifdef NDEBUG
        readOption(debug, "veryVerbose", debugParameters::veryVerbose);
    #else
        debugParameters::veryVerbose = true;
    #endif
    debugSourcePoint = -1;
    readOption(debug, "sourcePoint", debugSourcePoint);
    readOption(debug, "sourceTime", debugSourceTime);

    // CVODE options
    const object& cvode = conf.attr("cvodeTolerances");
    readOption(cvode, "relativeTolerance", integratorRelTol);
    readOption(cvode, "momentumAbsTol", integratorMomentumAbsTol);
    readOption(cvode, "energyAbsTol", integratorEnergyAbsTol);
    readOption(cvode, "speciesAbsTol", integratorSpeciesAbsTol);
    readOption(cvode, "minimumTimestep", integratorMinTimestep);

    // QSS integrator options
    const object& qss = conf.attr("qssTolerances");
    readOption(qss, "epsmax", qss_epsmax);
    readOption(qss, "epsmin", qss_epsmin);
    readOption(qss, "dtmin", qss_dtmin);
    readOption(qss, "dtmax", qss_dtmax);
    readOption(qss, "iterationCount", qss_iterationCount);
    readOption(qss, "abstol", qss_abstol);
    readOption(qss, "minval", qss_minval);
    readOption(qss, "stabilityCheck", qss_stabilityCheck);

    // Output file options
    const object& output = conf.attr("outputFiles");
    readOption(output, "saveProfiles", outputProfiles);
    if (outputProfiles) {
        readOption(output, "heatReleaseRate", outputHeatReleaseRate);
        readOption(output, "auxiliaryVariables", outputAuxiliaryVariables);
        readOption(output, "timeDerivatives", outputTimeDerivatives);
        readOption(output, "extraVariables", outputExtraVariables);
        outputFileNumber = 0;
        fileNumberOverride = readOption(output, "firstFileNumber", outputFileNumber);
    }
    readOption(output, "debugIntegratorStages", outputDebugIntegratorStages);

    // Termination Conditions
    const object& termination = conf.attr("terminationCondition");
    readOption(termination, "tEnd", tEnd);

    string terminationMeasurement;
    readOption(termination, "measurement", terminationMeasurement);
    terminateForSteadyQdot = (terminationMeasurement == "Q");

    if (terminateForSteadyQdot) {
        readOption(termination, "steadyPeriod", terminationPeriod);
        readOption(termination, "tolerance", terminationTolerance);
        readOption(termination, "abstol", terminationAbsTol);
    }

    // Create output directory if necessary
    if (!boost::filesystem::exists(outputDir)) {
        boost::filesystem::create_directory(outputDir);
    }

    if (boost::filesystem::exists(inputDir + "/" + gasMechanismFile)) {
        gasMechanismFile = inputDir + "/" + gasMechanismFile;
    }

    gridAlpha = (curvedFlame) ? 1 : 0;
}
