#pragma once
#include <string>
#include <iostream>

class debugParameters
{
public:
	static bool debugAdapt;
	static bool debugRegrid;
	static bool debugSundials;
	static bool debugJacobian;
	static bool debugCalcIC;
	static bool debugTimesteps;
	static bool debugSolverStats;
};

namespace debugType {
	enum debugType {
		adaptation,
		regridding,
		sundials
	};
};

void debugWrite(debugType::debugType type, std::string message);
