#pragma once

class configOptions
{
public:
	std::string inputDir;
	std::string outputDir;
	std::string restartFile;

	bool curvedDomain;

	bool overrideTu;
	bool overrideReactants;
	bool haveRestartFile;

	bool fixedBurnedVal;
	bool fixedLeftLoc;

	int regridStepInterval;
	int outputStepInterval;
	double regridTimeInterval;
	double outputTimeInterval;
	double maxTimestep;
	
	double idaRelTol;
	double idaContinuityAbsTol;
	double idaMomentumAbsTol;
	double idaEnergyAbsTol;
	double idaSpeciesAbsTol;
};
