#include "debugUtils.h"

bool debugParameters::debugAdapt = true;
bool debugParameters::debugRegrid = true;
bool debugParameters::debugSundials = false;


void debugWrite(debugType::debugType type, std::string message)
{
	if ((type == debugType::adaptation && debugParameters::debugAdapt) ||
		(type == debugType::regridding && debugParameters::debugRegrid) ||
		(type == debugType::sundials && debugParameters::debugSundials))
	{
		std::cout << message << std::endl;
	}
}