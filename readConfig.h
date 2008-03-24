#pragma once

class configOptions
{
public:
	std::string inputDir;
	std::string outputDir;
	std::string restartFile;

	bool overrideTu;
	bool overrideReactants;
	bool haveRestartFile;

	bool fixedBurnedVal;
	bool fixedLeftLoc;
	bool twinFlame;
	bool curvedFlame;

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


