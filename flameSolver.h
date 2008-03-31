#pragma once

#include "flameSys.h"

class flameSolver
{
public:
	void setOptions(const configOptions& options);
	void run(void);

	void calculateReactantMixture(void);
	bool checkTerminationCondition(void);

	// Time-series data
	dvector timeVector; // [s]
	dvector timestepVector; // [s]
	dvector heatReleaseRate; // [W/m^2]
	dvector consumptionSpeed; // [m/s]
	dvector flamePosition; // [m]

	configOptions options;
	flameSys theSys;
private:
    
};
