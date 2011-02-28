#include "readConfig.h"
#include "boost/filesystem.hpp"
#include "debugUtils.h"

#include <boost/python.hpp>

void configOptions::readOptionsFile(const std::string& filename)
{
    theConfig = new libconfig::Config;
    libconfig::Config& cfg = *theConfig;

    // read input file
    if (boost::filesystem::exists(filename)) {
        cfg.readFile(filename.c_str());
        logFile.write(format("Reading configuration options from %s") % filename);
    } else {
        throw debugException("configOptions::readOptionsFile: Error: Input file \"" + filename + "\" does not exist.");
    }


    cfg.setAutoConvert(true);

    // Read options from the configuration file

    // Paths
    readOption("paths.inputDir", inputDir, "input");
    readOption("paths.outputDir", outputDir, "output");

    // Chemistry
    bool haveMech = readOption("chemistry.mechanismFile", gasMechanismFile, "gri30.xml");
    if (haveMech) {
        readOption("chemistry.phaseID", gasPhaseID, "gas");
    } else {
        readOption("chemistry.phaseID", gasPhaseID, "gri30_multi");
    }

    std::string transportModel;
    readOption("chemistry.transportModel", transportModel, "Multi");
    if (transportModel=="Multi") {
        usingMultiTransport = true;
    } else if (transportModel=="Mix") {
        usingMultiTransport = false;
    } else {
        throw debugException("configOptions::readOptionsFile: Invalid Transport Model specified (general.transportModel).");
    }

    // Chemkin/AdapChem
    usingAdapChem = readOption("adapchem.mechanism", chemkinMechanismFile, "");
    if (usingAdapChem) {
        readOption("adapchem.input", adapchemInputFile, "adapchem.in");
        readOption("adapchem.models", adapchemModelsFile, "models.chem");
        readOption("adapchem.defaultModel", adapchemDefaultModelFile, "default.model");
        readOption("adapchem.donemodels", adapchemDonemodelsFile, "donemodels");
        readOption("adapchem.restart", adapchemRestartFile, "full.rstrt");
    }
    // Transport elimination
    readOption("adapchem.atol", adapchem_atol, 1e-8);
    readOption("adapchem.transportEliminationDiffusion", transportEliminationDiffusion, false);
    readOption("adapchem.transportEliminationConvection", transportEliminationConvection, false);
    readOption("adapchem.transportEliminationStepInterval", transportEliminationStepInterval, 1);

    // Grid
    readOption("initialCondition.nPoints", nPoints, 50);
    readOption("initialCondition.xLeft", xLeft, -0.005);
    readOption("initialCondition.xRight", xRight, 0.005);

    // Initial Condition
    haveRestartFile = readOptionQuietDefault("initialCondition.file", restartFile, "");
    useRelativeRestartPath = true;
    if (haveRestartFile) {
        haveRestartFile = boost::filesystem::exists(inputDir + "/" + restartFile);
        if (!haveRestartFile) {
            logFile.write(format("  Warning: couldn't find restart file '%s/%s'") %
                    inputDir % restartFile);
        }
    }

    overrideTu = readOption("initialCondition.Tu", Tu, 300);
    overrideReactants = readOption("initialCondition.fuel", fuel, "CH4:1.0");
    readOption("initialCondition.oxidizer", oxidizer, "O2:1.0, N2:3.76");
    readOption("initialCondition.equivalenceRatio", equivalenceRatio, 0.6);
    readOption("initialCondition.pressure", pressure, Cantera::OneAtm);


    if (readOptionQuietDefault("positionControl.xStag", xStag, 0)) {
        xStagControl = true;
    } else {
        xStagControl = false;
    }

    if (readOptionQuietDefault("positionControl.xInitial", xFlameInitial, 0.005) &&
            readOptionQuietDefault("positionControl.xFinal", xFlameFinal, 0.005) &&
            readOptionQuietDefault("positionControl.tStart", xFlameT0, 0.0) &&
            readOptionQuietDefault("positionControl.dt", xFlameDt, 1e-3)) {

        readOption("positionControl.integralGain",xFlameIntegralGain, 500);
        readOption("positionControl.proportionalGain",xFlameProportionalGain, 100);
        xFlameControl = true;
        xStagControl = true;
    } else {
        xFlameControl = false;
    }

    readOption("grid.centerGridMin", centerGridMin, 1.0e-4);

    readOption("grid.adaptation.vtol", vtol, 0.06);
    readOption("grid.adaptation.dvtol", dvtol, 0.2);

    readOption("grid.adaptation.rmTol", rmTol, 0.67);
    readOption("grid.adaptation.dampConst", dampConst, 0.5);
    readOption("grid.adaptation.gridMin", gridMin, 1.0e-6);
    readOption("grid.adaptation.gridMax", gridMax, 1.0e-3);
    readOption("grid.adaptation.uniformityTol", uniformityTol, 3.0);
    readOption("grid.adaptation.absvtol", absvtol, 1.0e-10);

    readOption("grid.regridding.boundaryTol", boundaryTol, 2.0e-5);
    readOption("grid.regridding.boundaryTolRm", boundaryTolRm, 5.0e-6);
    readOption("grid.regridding.addPointCount", addPointCount, 2);

    haveTStart = readOption("times.tStart", tStart, 0.0e0);
    readOption("terminationCondition.tEnd", tEnd, 10);

    readOption("general.fixedBurnedVal", fixedBurnedVal, true);
    readOption("general.fixedLeftLocation", fixedLeftLoc, false);
    readOption("general.unburnedLeft",unburnedLeft, true);
    readOption("general.curvedFlame",curvedFlame, false);
    readOption("general.twinFlame",twinFlame, false);

    if (cfg.exists("times.regridTimeInterval") || cfg.exists("times.regridStepInterval")) {
        readOptionQuietDefault("times.regridTimeInterval",regridTimeInterval, 100);
        readOptionQuietDefault("times.regridStepInterval",regridStepInterval, 100000);
    } else {
        readOption("times.regridTimeInterval",regridTimeInterval, 0.005);
        readOption("times.regridStepInterval",regridStepInterval, 100);
    }

    if (cfg.exists("times.outputTimeInterval") || cfg.exists("times.outputStepInterval")) {
        readOptionQuietDefault("times.outputTimeInterval", outputTimeInterval, 100);
        readOptionQuietDefault("times.outputStepInterval", outputStepInterval, 100000);
    } else {
        readOption("times.outputTimeInterval", outputTimeInterval, 0.001);
        readOption("times.outputStepInterval", outputStepInterval, 10);
    }

    if (cfg.exists("times.profileTimeInterval") || cfg.exists("times.profileStepInterval")) {
        readOptionQuietDefault("times.profileTimeInterval", profileTimeInterval, 100);
        readOptionQuietDefault("times.profileStepInterval", profileStepInterval, 100000);
    } else {
        readOption("times.profileTimeInterval", profileTimeInterval, 0.050);
        readOption("times.profileStepInterval", profileStepInterval, 2000);
    }

    readOption("times.currentStateStepInterval", currentStateStepInterval, 2000);
    readOption("times.terminateStepInterval", terminateStepInterval, 100);
    readOption("times.maxTimestep",maxTimestep, 1.0e-3);
    readOption("times.globalTimestep", globalTimestep, 1.0e-4);
    readOption("times.diffusionTimestep", diffusionTimestep, 1.0e-6);

    readOption("debug.adaptation",debugParameters::debugAdapt, false);
    readOption("debug.regridding",debugParameters::debugRegrid, false);
    readOption("debug.sundials",debugParameters::debugSundials, false);
    readOption("debug.jacobian",debugParameters::debugJacobian, false);
    readOption("debug.calcIC",debugParameters::debugCalcIC, false);
    readOption("debug.timesteps",debugParameters::debugTimesteps, true);
    readOption("debug.solverStats",debugParameters::debugSolverStats, true);
    readOption("debug.performanceStats",debugParameters::debugPerformanceStats, true);
    readOption("debug.flameRadiusControl",debugParameters::debugFlameRadiusControl, false);
    #ifdef NDEBUG
        readOption("debug.veryVerbose", debugParameters::veryVerbose, false);
    #else
        debugParameters::veryVerbose = true;
    #endif

    readOption("integrator.chemistryIntegrator", chemistryIntegrator, "cvode");
    if (chemistryIntegrator != "cvode" && chemistryIntegrator != "qss") {
        throw debugException("Unknown chemistryIntegrator: '" + chemistryIntegrator + "'.");
    }
    readOption("integrator.relativeTolerance", integratorRelTol, 1e-5);
    readOption("integrator.momentumAbsTol", integratorMomentumAbsTol, 1e-8);
    readOption("integrator.energyAbsTol", integratorEnergyAbsTol, 1e-6);
    readOption("integrator.speciesAbsTol", integratorSpeciesAbsTol, 1e-10);
    readOption("integrator.minimumTimestep", integratorMinTimestep, 1e-18);

    readOption("integrator.qss.epsmax", qss_epsmax, 10);
    readOption("integrator.qss.epsmin", qss_epsmin, 0.02);
    readOption("integrator.qss.dtmin", qss_dtmin, 1e-16);
    readOption("integrator.qss.iterationCount", qss_iterationCount, 1);
    readOption("integrator.qss.abstol", qss_abstol, 1e-11);
    readOption("integrator.qss.minval", qss_minval, 1e-30);
    readOption("integrator.qss.stabilityCheck", qss_stabilityCheck, false);

    if (readOption("outputFiles.saveProfiles", outputProfiles, true)) {
        readOption("outputFiles.heatReleaseRate", outputHeatReleaseRate, true);
        readOption("outputFiles.auxiliaryVariables", outputAuxiliaryVariables, false);
        readOption("outputFiles.timeDerivatives", outputTimeDerivatives, false);
        readOption("outputFiles.residualComponents", outputResidualComponents, false);
        fileNumberOverride = readOption("outputFiles.firstFileNumber", outputFileNumber, 0);
    }
    readOption("outputFiles.splitHeatReleaseRate", outputSplitHeatReleaseRate, false);
    readOption("outputFiles.debugIntegratorStages", outputDebugIntegratorStages, false);
    stopIfError = readOption("general.errorStopCount",errorStopCount, 1);

    readOptionQuietDefault("debug.sourcePoint", debugSourcePoint, -1);
    readOptionQuietDefault("debug.debugSourceTime", debugSourceTime, 0);

    if (cfg.exists("strainParameters.list")) {
        libconfig::Setting& strainSetting = cfg.lookup("strainParameters.list");
        int strainCount = strainSetting.getLength();
        for (int i=0; i<strainCount; i++) {
            strainRateList.push_back(strainSetting[i]);
        }
        logFile.write("    read option: strainParameters.list = ", false);
        logFile.write(strainRateList);
        multiRun = true;
    } else {
        // Strain Rate Parameters
        readOption("strainParameters.initial", strainRateInitial, 100);
        readOption("strainParameters.final", strainRateFinal, 100);
        readOption("strainParameters.tStart", strainRateT0, 1.0e-3);
        readOption("strainParameters.dt", strainRateDt, 0.0e0);
        multiRun = false;
    }

    if (multiRun) {
        terminateForSteadyQdot = true;
    } else {
        std::string terminationMeasurement;
        readOptionQuietDefault("terminationCondition.measurement",terminationMeasurement, "");
        terminateForSteadyQdot = (terminationMeasurement == "Q");
    }

    if (terminateForSteadyQdot) {
        readOption("terminationCondition.time",terminationPeriod, 0.01);
        readOption("terminationCondition.timeMax",terminationMaxTime, 2.0);
        readOption("terminationCondition.tolerance", terminationTolerance, 1e-4);
        readOption("terminationCondition.abstol",terminationAbsTol, 0.5);
    }

    if (!boost::filesystem::exists(outputDir)) {
        boost::filesystem::create_directory(outputDir);
    }

    if (boost::filesystem::exists(inputDir + "/" + gasMechanismFile)) {
        gasMechanismFile = inputDir + "/" + gasMechanismFile;
    }

    gridAlpha = (curvedFlame) ? 1 : 0;

    logFile.write("Finished reading configuration options.");
    delete theConfig;
}

template <class T1, class T2>
bool configOptions::readOption(const std::string name, T1& value, const T2 defaultVal)
{
    bool readVal = theConfig->lookupValue(name, value);
    if (readVal) {
        logFile.write(format("    read option: %s = %s") % name % value);
    } else {
        value = defaultVal;
        logFile.write(format(" * used default: %s = %s") % name % defaultVal);
    }
    return readVal;
}

template <class T1, class T2>
bool configOptions::readOptionQuietDefault(const std::string name, T1& value, const T2 defaultVal)
{
    bool readVal = theConfig->lookupValue(name, value);
    if (readVal) {
        logFile.write(format("    read option: %s = %s") % name % value);
    } else {
        value = defaultVal;
    }
    return readVal;
}

const size_t kMomentum = 0;
const size_t kEnergy = 1;
const size_t kSpecies = 2;
const size_t kWmx = 2; // never used in the same systems as kSpecies

template <class T1>
bool configOptions::readOption(const boost::python::object& conf, const char* name, T1& value)
{
    const boost::python::object& obj = conf.attr(name);
    if (obj.ptr() != Py_None) {
        value = boost::python::extract<T1>(obj);
        return true;
    } else {
        return false;
    }
}

configOptions::configOptions(const boost::python::object& conf)
{
    using boost::python::extract;
    using boost::python::object;
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

    readOption(general, "curvedFlame", curvedFlame);
    readOption(general, "twinFlame", twinFlame);
    stopIfError = readOption(general, "errorStopCount", errorStopCount);

    // Chemistry
    const object& chem = conf.attr("chemistry");
    readOption(chem, "mechanismFile", gasMechanismFile);
    readOption(chem, "phaseID", gasPhaseID);
    string transportModel = extract<string>(chem.attr("transportModel"));
    if (transportModel == "Multi") {
        usingMultiTransport = true;
    } else if (transportModel == "Mix") {
        usingMultiTransport = false;
    } else {
        throw debugException("configOptions::readOptionsFile: Invalid Transport Model specified (general.transportModel).");
    }

    // Chemkin//Adapchem
    const object& adapchem = conf.attr("adapchem");
    if (adapchem != None) {
        usingAdapChem = true;
        readOption(adapchem, "inputFile", adapchemInputFile);
        readOption(adapchem, "models", adapchemModelsFile);
        readOption(adapchem, "defaultModel", adapchemDefaultModelFile);
        readOption(adapchem, "donemodels", adapchemDefaultModelFile);
        readOption(adapchem, "restart", adapchemRestartFile);
    } else {
        usingAdapChem = false;
    }

    // Transport species elimination
    const object& transport = conf.attr("transportElimination");
    readOption(transport, "atol", adapchem_atol);
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
    if (haveRestartFile) {
        haveRestartFile = boost::filesystem::exists(inputDir + "/" + restartFile);
        if (!haveRestartFile) {
            logFile.write(format("WARNING: couldn't find restart file '%s/%s'") %
                inputDir % restartFile);
        }
    }

    overrideTu = readOption(ic, "Tu", Tu);
    overrideReactants = readOption(ic, "fuel", fuel);
    readOption(ic, "flameType", flameType);
    readOption(ic, "oxidizer", oxidizer);
    readOption(ic, "Tfuel", Tfuel);
    readOption(ic, "Toxidizer", Toxidizer);
    readOption(ic, "equivalenceRatio", equivalenceRatio);
    readOption(ic, "pressure", pressure);

    // Strain Rate Parameters
    const object& strain = conf.attr("strainParameters");
    if (strain.attr("rates") != None) {
        readOption(strain, "rates", strainRateList);
        multiRun = true;
    } else {
        readOption(strain, "initial", strainRateInitial);
        readOption(strain, "final", strainRateFinal);
        readOption(strain, "tStart", strainRateT0);
        readOption(strain, "dt", strainRateDt);
        multiRun = false;
    }

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
    readOption(times, "maxTimestep", maxTimestep);
    readOption(times, "globalTimestep", globalTimestep);
    readOption(times, "diffusionTimestep", diffusionTimestep);

    // Debug parameters
    const object& debug = conf.attr("debug");
    readOption(debug, "adaptation", debugParameters::debugAdapt);
    readOption(debug, "regridding", debugParameters::debugRegrid);
    readOption(debug, "sundials", debugParameters::debugSundials);
    readOption(debug, "jacobian", debugParameters::debugJacobian);
    readOption(debug, "calcIC", debugParameters::debugCalcIC);
    readOption(debug, "timesteps", debugParameters::debugTimesteps);
    readOption(debug, "solverStats", debugParameters::debugSolverStats);
    readOption(debug, "performanceStats", debugParameters::debugPerformanceStats);
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
    if (cvode != None) {
        readOption(cvode, "relativeTolerance", integratorRelTol);
        readOption(cvode, "momentumAbsTol", integratorMomentumAbsTol);
        readOption(cvode, "energyAbsTol", integratorEnergyAbsTol);
        readOption(cvode, "speciesAbsTol", integratorSpeciesAbsTol);
        readOption(cvode, "minimumTimestep", integratorMinTimestep);
    }

    // QSS integrator options
    const object& qss = conf.attr("qssTolerances");
    readOption(qss, "epsmax", qss_epsmax);
    readOption(qss, "epsmin", qss_epsmin);
    readOption(qss, "dtmin", qss_dtmin);
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
        readOption(output, "residualComponents", outputResidualComponents);
        outputFileNumber = 0;
        fileNumberOverride = readOption(output, "firstFileNumber", outputFileNumber);
    }
    readOption(output, "splitHeatReleaseRate", outputSplitHeatReleaseRate);
    readOption(output, "debugIntegratorStages", outputDebugIntegratorStages);

    // Termination Conditions
    const object& termination = conf.attr("terminationCondition");
    readOption(termination, "tEnd", tEnd);
    if (multiRun) {
        terminateForSteadyQdot = true;
    } else {
        string terminationMeasurement;
        readOption(termination, "measurement", terminationMeasurement);
        terminateForSteadyQdot = (terminationMeasurement == "Q");
    }

    if (terminateForSteadyQdot) {
        readOption(termination, "steadyPeriod", terminationPeriod);
        readOption(termination, "timeMax", terminationMaxTime);
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
