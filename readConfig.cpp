#include "strainedFlameSys.h"
#include "libconfig.h++"
#include "boost/filesystem.hpp"
#include "debugUtils.h"

void configOptions::readOptionsFile(const std::string& filename)
{
	// read input file
	libconfig::Config cfg;
	if (boost::filesystem::exists(filename)) {
		cfg.readFile(filename.c_str());
		cout << "Reading configuration options from " << filename << endl;
	} else {
		cout << "readOptionsFile: Error: Input file \"" << filename << "\" does not exist." << endl;
		throw;
	}

	cfg.setAutoConvert(true);

	// These are default values for the configuration options:

	// Paths
	inputDir = "input";
	outputDir = "output";

	// Chemistry
	gasMechanismFile = "gri30.xml";
	gasPhaseID = "gri30_multi";

	// Initial Conditions
	restartFile = "";
	reactants = "O2:1.0, N2:3.76, CH4:0.5";
	diluent = "Ar:1.0";
	Tu = 300;
	Tb = 2000;
	pressure = Cantera::OneAtm;

	// Strain Rate Parameters
	strainRateInitial = 100;
	strainRateFinal = 100;
	strainRateDt = 1e-3;
	strainRateT0 = 0;
	
	curvedFlame = false;
	fixedLeftLoc = false;
	twinFlame = false;

	// Grid	
	nPoints = 50;
	xLeft = -0.05;
	xRight = 0.05;

	vtol = 0.04;
	dvtol = 0.4;
	rmTol = 0.67;
	dampConst = 5000;
	gridMin = 2.0e-4;
	gridMax = 0.2;
	uniformityTol = 2.7;
	absvtol = 1e-10;
	
	boundaryTol = 2e-5;
	boundaryTolRm = 5e-6;

	fixedBurnedVal = true;
	unburnedLeft = true;
	addPointCount = 1;

	// Times
	tStart = 0;
	tEnd = 0.05;
	maxTimestep = 10000;

	regridTimeInterval = 100;
	outputTimeInterval = 100;
	regridStepInterval = 10000000;
	outputStepInterval = 10000000;

	idaRelTol = 1e-5;
	idaContinuityAbsTol = 1e-6;
	idaMomentumAbsTol = 1e-6;
	idaEnergyAbsTol = 1e-6;
	idaSpeciesAbsTol = 1e-10;
	
	// Flags


	// Read options from the configuration file
	cfg.lookupValue("paths.inputDir",inputDir);
	cfg.lookupValue("paths.outputDir",outputDir);
	
	cfg.lookupValue("chemistry.mechanismFile",gasMechanismFile);
	cfg.lookupValue("chemistry.phaseID",gasPhaseID);

	cfg.lookupValue("grid.nPoints",nPoints);
	cfg.lookupValue("grid.xLeft",xLeft);
	cfg.lookupValue("grid.xRight",xRight);

	haveRestartFile = cfg.lookupValue("InitialCondition.file",restartFile);
	overrideTu = cfg.lookupValue("InitialCondition.Tu",Tu);
	overrideReactants = cfg.lookupValue("InitialCondition.reactants",reactants);
	cfg.lookupValue("InitialCondition.diluent",diluent);
	cfg.lookupValue("InitialCondition.pressure",pressure);

	cfg.lookupValue("StrainParameters.initial",strainRateInitial);
	cfg.lookupValue("StrainParameters.final",strainRateFinal);
	cfg.lookupValue("StrainParameters.start",strainRateT0);
	cfg.lookupValue("StrainParameters.dt",strainRateDt);

	cfg.lookupValue("grid.adaptation.vtol",vtol);
	cfg.lookupValue("grid.adaptation.dvtol",dvtol);
	cfg.lookupValue("grid.adaptation.rmTol",rmTol);
	cfg.lookupValue("grid.adaptation.dampConst",dampConst);
	cfg.lookupValue("grid.adaptation.gridMin",gridMin);
	cfg.lookupValue("grid.adaptation.gridMax",gridMax);
	cfg.lookupValue("grid.adaptation.uniformityTol",uniformityTol);
	cfg.lookupValue("grid.adaptation.absvtol",absvtol);
	
	cfg.lookupValue("grid.regridding.boundaryTol",boundaryTol);
	cfg.lookupValue("grid.regridding.boundaryTolRm",boundaryTolRm);
	cfg.lookupValue("grid.regridding.addPointCount",addPointCount);

	cfg.lookupValue("times.tStart",tStart);
	cfg.lookupValue("times.tEnd",tEnd);

	cfg.lookupValue("general.fixedBurnedVal",fixedBurnedVal);
	cfg.lookupValue("general.fixedLeftLocation",fixedLeftLoc);
	cfg.lookupValue("general.unburnedLeft",unburnedLeft);
	cfg.lookupValue("general.curvedFlame",curvedFlame);
	cfg.lookupValue("general.twinFlame",twinFlame);

	cfg.lookupValue("times.regridTimeInterval",regridTimeInterval);
	cfg.lookupValue("times.regridStepInterval",regridStepInterval);
	cfg.lookupValue("times.outputTimeInterval",outputTimeInterval);
	cfg.lookupValue("times.outputStepInterval",outputStepInterval);
	cfg.lookupValue("times.maxTimestep",maxTimestep);

	cfg.lookupValue("debug.adaptation",debugParameters::debugAdapt);
	cfg.lookupValue("debug.regridding",debugParameters::debugRegrid);
	cfg.lookupValue("debug.sundials",debugParameters::debugSundials);
	cfg.lookupValue("debug.jacobian",debugParameters::debugJacobian);
	cfg.lookupValue("debug.calcIC",debugParameters::debugCalcIC);

	cfg.lookupValue("integrator.relativeTolerance",idaRelTol);
	cfg.lookupValue("integrator.continuityAbsTol",idaContinuityAbsTol);
	cfg.lookupValue("integrator.momentumAbsTol",idaMomentumAbsTol);
	cfg.lookupValue("integrator.energyAbsTol",idaEnergyAbsTol);
	cfg.lookupValue("integrator.speciesAbsTol",idaSpeciesAbsTol);

	if (haveRestartFile) {
		haveRestartFile = boost::filesystem::exists(inputDir + "/" + restartFile);
	}

	if (!boost::filesystem::exists(outputDir)) {
		boost::filesystem::create_directory(outputDir);
	}

	if (boost::filesystem::exists(inputDir + "/" + gasMechanismFile)) {
		gasMechanismFile = inputDir + "/" + gasMechanismFile;
	}

	gridAlpha = (curvedFlame) ? 1 : 0;

	kContinuity = 0;
	kMomentum = 1;
	kEnergy = 2;
	kSpecies = 3;

	// If neither step nor time intervals have been specified, use a default step interval
	if (outputTimeInterval == 100 && outputStepInterval == 10000000) {
		outputStepInterval = 50;
	}

	if (regridTimeInterval == 100 && regridStepInterval == 10000000) {
		regridStepInterval = 20;
	}
}
