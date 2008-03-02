#include "strainedFlameSys.h"
#include "libconfig.h++"
#include "boost/filesystem.hpp"

void strainedFlameSys::readOptionsFile(std::string filename)
{
	// read input file
	libconfig::Config configFile;
	if (boost::filesystem::exists(filename)) {
		configFile.readFile(filename.c_str());
		cout << "Reading configuration opions from " << filename << endl;
	} else {
		cout << "readOptionsFile: Error: Input file \"" << filename << "\" does not exist." << endl;
		throw;
	}

	configFile.setAutoConvert(true);

	// These are default values for the configuration options

	// Paths
	options.inputDir = "input";
	options.outputDir = "output";

	// Chemistry
	gas.mechanismFile = "gri30.xml";
	gas.phaseID = "gri30_multi";
	reactants = "O2:1.0, N2:3.76, CH4:0.5";
	diluent = "Ar:1.0";

	// Flow
	Tu = 300;
	Tb = 2000;
	strainRate = 100;
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
	configFile.lookupValue("paths.inputDir",options.inputDir);
	configFile.lookupValue("paths.outputDir",options.outputDir);

	configFile.lookupValue("chemistry.mechanismFile",gas.mechanismFile);
	configFile.lookupValue("chemistry.phaseID",gas.phaseID);
	configFile.lookupValue("chemistry.reactants",reactants);
	configFile.lookupValue("chemistry.diluent",diluent);
	configFile.lookupValue("chemistry.pressure",gas.pressure);

	configFile.lookupValue("grid.nPoints",nPoints);
	configFile.lookupValue("grid.xLeft",xLeft);
	configFile.lookupValue("grid.xRight",xRight);

	configFile.lookupValue("flow.Tu",Tu);
	configFile.lookupValue("flow.Tb",Tb);

	configFile.lookupValue("flow.rhou",rhou);
	configFile.lookupValue("flow.strainRate",strainRate);

	configFile.lookupValue("tStart",tStart);
	configFile.lookupValue("tEnd",tEnd);

	configFile.lookupValue("adaptation.vtol",grid.vtol);
	configFile.lookupValue("adaptation.dvtol",grid.dvtol);
	configFile.lookupValue("adaptation.rmTol",grid.rmTol);
	configFile.lookupValue("adaptation.dampConst",grid.dampConst);
	configFile.lookupValue("adaptation.gridMin",grid.gridMin);
	configFile.lookupValue("adaptation.gridMax",grid.gridMax);
	configFile.lookupValue("adaptation.uniformityTol",grid.uniformityTol);

	configFile.lookupValue("general.fixedBurnedValFlag",grid.fixedBurnedValFlag);
	configFile.lookupValue("general.fixedLeftLocation",grid.fixedLeftLoc);
	configFile.lookupValue("general.unburnedLeft",grid.unburnedLeft);

	configFile.lookupValue("regridding.boundaryTol",grid.boundaryTol);
	configFile.lookupValue("regridding.boundaryTolRm",grid.boundaryTolRm);


	grid.alpha = (options.curvedDomain) ? 1 : 0;

	grid.kContinuity = 0;
	grid.kMomentum = 1;
}