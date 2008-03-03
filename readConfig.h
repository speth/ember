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
};
