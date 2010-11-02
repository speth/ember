#include "readConfig.h"
#include "boost/filesystem.hpp"
#include "debugUtils.h"

using std::cout;
using std::endl;

void configOptions::readOptionsFile(const std::string& filename)
{
    theConfig = new libconfig::Config;
    libconfig::Config& cfg = *theConfig;

    // read input file
    if (boost::filesystem::exists(filename)) {
        cfg.readFile(filename.c_str());
        cout << "Reading configuration options from " << filename << endl;
    } else {
        throw debugException("configOptions::readOptionsFile: Error: Input file \"" + filename + "\" does not exist.");
    }

    cout << std::boolalpha; // prints "true" and "false" rather than 1 and 0
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
        readOption("adapchem.atol", adapchem_atol, 1e-8);
    }
    readOption("adapchem.transportElimination", transportElimination, false);

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
            cout << "  Warning: couldn't find restart file \"" << inputDir+"/"+restartFile << "\"" << endl;
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
    readOption("grid.adaptation.vtolCont", vtolCont, vtol);
    readOption("grid.adaptation.dvtolCont", dvtolCont, dvtol);

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
    readOption("general.centeredDifferences",centeredDifferences, false);
    readOption("general.steadyOnly",steadyOnly, false);

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
    readOption("times.integratorRestartInterval",integratorRestartInterval, 200);
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

    readOption("integrator.relativeTolerance", idaRelTol, 1e-5);
    readOption("integrator.relativeToleranceLow", idaRelTolLow, 1e-3);
    readOption("integrator.continuityAbsTol", idaContinuityAbsTol, 1e-8);
    readOption("integrator.momentumAbsTol", idaMomentumAbsTol, 1e-8);
    readOption("integrator.energyAbsTol", idaEnergyAbsTol, 1e-6);
    readOption("integrator.speciesAbsTol", idaSpeciesAbsTol, 1e-10);

    if (readOption("outputFiles.saveProfiles", outputProfiles, true)) {
        readOption("outputFiles.heatReleaseRate", outputHeatReleaseRate, true);
        readOption("outputFiles.auxiliaryVariables", outputAuxiliaryVariables, false);
        readOption("outputFiles.timeDerivatives", outputTimeDerivatives, false);
        readOption("outputFiles.residualComponents", outputResidualComponents, false);
        fileNumberOverride = readOption("outputFiles.firstFileNumber", outputFileNumber, 0);
    }
    stopIfError = readOption("general.errorStopCount",errorStopCount, 1);

    if (cfg.exists("strainParameters.list")) {
        libconfig::Setting& strainSetting = cfg.lookup("strainParameters.list");
        int strainCount = strainSetting.getLength();
        for (int i=0; i<strainCount; i++) {
            strainRateList.push_back(strainSetting[i]);
        }
        cout << "    read option: strainParameters.list = " << strainRateList << endl;
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
        readOption("terminationCondition.timeMax",terminationMaxTime, 2.0);
        readOption("terminationCondition.timeLow",terminationPeriodLow, 0.04);
        readOption("terminationCondition.timeHigh",terminationPeriodHigh, 0.10);
        readOption("terminationCondition.tolerance", terminationTolerance, 1e-4);
        readOption("terminationCondition.toleranceLow", terminationToleranceLow, 1e-4);
        readOption("terminationCondition.abstol",terminationAbsTol, 0.5);
    } else {
        std::string terminationMeasurement;
        readOptionQuietDefault("terminationCondition.measurement",terminationMeasurement, "");
        terminateForSteadyQdot = (terminationMeasurement == "Q");
        if (terminateForSteadyQdot) {
            readOption("terminationCondition.timeMax",terminationMaxTime, 2.0);
            readOption("terminationCondition.time",terminationPeriod, 0.10);
            readOption("terminationCondition.tolerance", terminationTolerance, 1e-4);
            readOption("terminationCondition.abstol",terminationAbsTol, 0.5);
        }
    }

    if (!boost::filesystem::exists(outputDir)) {
        boost::filesystem::create_directory(outputDir);
    }

    if (boost::filesystem::exists(inputDir + "/" + gasMechanismFile)) {
        gasMechanismFile = inputDir + "/" + gasMechanismFile;
    }

    gridAlpha = (curvedFlame) ? 1 : 0;

    cout << "Finished reading configuration options." << endl;
    delete theConfig;
}

template <class T1, class T2>
bool configOptions::readOption(const std::string name, T1& value, const T2 defaultVal)
{
    bool readVal = theConfig->lookupValue(name, value);
    if (readVal) {
        cout << "    read option: " << name << " = " << value << endl;
    } else {
        value = defaultVal;
        cout << " * used default: " << name << " = " << defaultVal << endl;
    }
    return readVal;
}

template <class T1, class T2>
bool configOptions::readOptionQuietDefault(const std::string name, T1& value, const T2 defaultVal)
{
    bool readVal = theConfig->lookupValue(name, value);
    if (readVal) {
            cout << "    read option: " << name << " = " << value << endl;
    } else {
        value = defaultVal;
    }
    return readVal;
}

const size_t kMomentum = 0;
const size_t kEnergy = 1;
const size_t kSpecies = 2;
