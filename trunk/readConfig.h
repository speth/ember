#pragma once
#include <string>
#include "mathUtils.h"

class configOptions
{
public:
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

	int regridStepInterval;
	int outputStepInterval;
	int profileStepInterval;
	int integratorRestartInterval;
	double regridTimeInterval;
	double outputTimeInterval;
	double profileTimeInterval;
	double maxTimestep;

	double idaRelTol;
	double idaRelTolLow;
	double idaContinuityAbsTol;
	double idaMomentumAbsTol;
	double idaEnergyAbsTol;
	double idaSpeciesAbsTol;
	bool enforceNonnegativeSpecies;

	std::string gasMechanismFile;
	std::string gasPhaseID;
	bool usingMultiTransport;
	bool singleCanteraObject;

	std::string fuel;
	std::string oxidizer;
	double equivalenceRatio;
	dvector reactants; // mole fractions
	dvector products;
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
	bool outputResidualComponents;
	bool outputProfiles;

	bool terminateForSteadyQdot;
	double terminationTolerance; // relative tolerance
	double terminationToleranceLow;
	double terminationAbsTol; // absolute tolerance
	double terminationPeriod;
	double terminationPeriodLow;
	double terminationPeriodHigh;
	double terminationMaxTime;

	int errorStopCount;
	bool stopIfError;

	dvector strainRateList;

	int outputFileNumber; // number of output files written
	bool fileNumberOverride; // true if outputFileNumbe was given in the input file

    double centerGridMin;

	bool flameRadiusControl;
	double rFlameInitial, rFlameFinal;
	double rFlameDt, rFlameT0;
	double rFlameIntegralGain;
	double rFlameProportionalGain;
	double rFlameDerivativeGain;
	double rFlameUpdateStepInterval, rFlameUpdateTimeInterval;

	bool stagnationRadiusControl;
	double rStag;

	int numberOfThreads;
};
