#include "strainedFlameSys.h"
#include "libconfig.h++"
#include "boost/filesystem.hpp"

void strainedFlameSys::readOptionsFile(std::string filename)
{
	// read input file
	libconfig::Config cfg;
	if (boost::filesystem::exists(filename)) {
		cfg.readFile(filename.c_str());
		cout << "Reading configuration opions from " << filename << endl;
	} else {
		cout << "readOptionsFile: Error: Input file \"" << filename << "\" does not exist." << endl;
		throw;
	}

	cfg.setAutoConvert(true);

	// These are default values for the configuration options:

	// Paths
	options.inputDir = "input";
	options.outputDir = "output";

	// Chemistry
	gas.mechanismFile = "gri30.xml";
	gas.phaseID = "gri30_multi";

	// Initial Conditions
	options.restartFile = "";
	reactants = "O2:1.0, N2:3.76, CH4:0.5";
	diluent = "Ar:1.0";
	Tu = 300;
	Tb = 2000;

	// Strain Rate Parameters
	strainRateInitial = 100;
	strainRateFinal = 100;
	strainRateDt = 1e-3;
	strainRateT0 = 0;
	
	options.curvedDomain = false;

	// Grid	
	nPoints = 50;
	xLeft = -0.05;
	xRight = 0.05;

	grid.vtol = 0.04;
	grid.dvtol = 0.4;
	grid.rmTol = 0.67;
	grid.dampConst = 5000;
	grid.gridMin = 2.0e-4;
	grid.gridMax = 0.2;
	grid.uniformityTol = 2.7;
	
	grid.boundaryTol = 2e-5;
	grid.boundaryTolRm = 5e-6;

	grid.fixedBurnedValFlag = true;
	grid.fixedLeftLoc = false;
	grid.unburnedLeft = true;

	// Integrator
	tStart = 0;
	tEnd = 0.05;

	// Flags


	// Read options from the configuration file
	cfg.lookupValue("paths.inputDir",options.inputDir);
	cfg.lookupValue("paths.outputDir",options.outputDir);
	

	cfg.lookupValue("chemistry.mechanismFile",gas.mechanismFile);
	cfg.lookupValue("chemistry.phaseID",gas.phaseID);


	cfg.lookupValue("grid.nPoints",nPoints);
	cfg.lookupValue("grid.xLeft",xLeft);
	cfg.lookupValue("grid.xRight",xRight);

	options.haveRestartFile = cfg.lookupValue("InitialCondition.file",options.restartFile);
	options.overrideTu = cfg.lookupValue("InitialCondition.Tu",Tu);
	options.overrideReactants = cfg.lookupValue("InitialCondition.reactants",reactants);
	cfg.lookupValue("InitialCondition.diluent",diluent);
	cfg.lookupValue("InitialCondition.pressure",gas.pressure);

	cfg.lookupValue("StrainParameters.initial",strainRateInitial);
	cfg.lookupValue("StrainParameters.final",strainRateFinal);
	cfg.lookupValue("StrainParameters.start",strainRateT0);
	cfg.lookupValue("StrainParameters.dt",strainRateDt);

	cfg.lookupValue("grid.adaptation.vtol",grid.vtol);
	cfg.lookupValue("grid.adaptation.dvtol",grid.dvtol);
	cfg.lookupValue("grid.adaptation.rmTol",grid.rmTol);
	cfg.lookupValue("grid.adaptation.dampConst",grid.dampConst);
	cfg.lookupValue("grid.adaptation.gridMin",grid.gridMin);
	cfg.lookupValue("grid.adaptation.gridMax",grid.gridMax);
	cfg.lookupValue("grid.adaptation.uniformityTol",grid.uniformityTol);
	
	cfg.lookupValue("grid.regridding.boundaryTol",grid.boundaryTol);
	cfg.lookupValue("grid.regridding.boundaryTolRm",grid.boundaryTolRm);

	cfg.lookupValue("tStart",tStart);
	cfg.lookupValue("tEnd",tEnd);

	cfg.lookupValue("general.fixedBurnedValFlag",grid.fixedBurnedValFlag);
	cfg.lookupValue("general.fixedLeftLocation",grid.fixedLeftLoc);
	cfg.lookupValue("general.unburnedLeft",grid.unburnedLeft);

	if (options.haveRestartFile) {
		options.haveRestartFile = 
			boost::filesystem::exists(options.inputDir + "/" + options.restartFile);
	}

	if (!boost::filesystem::exists(options.outputDir)) {
		boost::filesystem::create_directory(options.outputDir);
	}

	grid.alpha = (options.curvedDomain) ? 1 : 0;

	grid.kContinuity = 0;
	grid.kMomentum = 1;
}