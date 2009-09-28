#include "debugUtils.h"

bool debugParameters::debugAdapt;
bool debugParameters::debugRegrid;
bool debugParameters::debugSundials;
bool debugParameters::debugJacobian;
bool debugParameters::debugCalcIC;
bool debugParameters::debugTimesteps;
bool debugParameters::debugSolverStats;
bool debugParameters::debugPerformanceStats;
bool debugParameters::debugFlameRadiusControl;

debugException::debugException(void)
{
    errorString = "debugException: unspecified error.";
}

debugException::debugException(const std::string error)
{
    errorString = error;
}
