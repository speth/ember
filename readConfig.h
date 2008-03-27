#pragma once

class configOptions
{
public:
	void readOptionsFile(const std::string& filename);

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

	std::string gasMechanismFile;
	std::string gasPhaseID;
	std::string reactants;
	std::string diluent;
	double pressure;

	double Tu, Tb;
	double strainRateInitial, strainRateFinal;
	double strainRateDt, strainRateT0;

	int nPoints;
	double xLeft, xRight;

	double vtol, dvtol, rmTol, dampConst, gridMin, gridMax;
	double uniformityTol, absvtol;
	double boundaryTol, boundaryTolRm;

	bool unburnedLeft;
	int addPointCount;

	double tStart, tEnd;

	int gridAlpha;
	int kContinuity, kMomentum, kEnergy, kSpecies;

	bool outputAuxiliaryVariables;
	bool outputTimeDerivatives;
	bool outputHeatReleaseRate;
};


