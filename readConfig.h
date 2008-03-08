#pragma once

class configOptions
{
public:
	string inputDir;
	string outputDir;
	string restartFile;

	bool curvedDomain;

	bool overrideTu;
	bool overrideReactants;
	bool haveRestartFile;

	int regridStepInterval;
	int outputStepInterval;
	double regridTimeInterval;
	double outputTimeInterval;
	
};
