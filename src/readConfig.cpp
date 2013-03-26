#include "readConfig.h"
#include "debugUtils.h"

const size_t kMomentum = 0;
const size_t kEnergy = 1;
const size_t kSpecies = 2;
const size_t kWmx = 2; // never used in the same systems as kSpecies

void ConfigOptions::setContinuityBC(const std::string& condition)
{
    if (condition == "fixedLeft") {
        continuityBC = ContinuityBoundaryCondition::Left;
    } else if (condition == "fixedRight") {
        continuityBC = ContinuityBoundaryCondition::Right;
    } else if (condition == "fixedQdot") {
        continuityBC = ContinuityBoundaryCondition::Qdot;
    } else if (condition == "fixedTemperature") {
        continuityBC = ContinuityBoundaryCondition::Temp;
    } else if (condition == "stagnationPoint") {
        continuityBC = ContinuityBoundaryCondition::Zero;
    }
}

bool ConfigOptions::debugIntegratorStages(double t) const
{
    return (outputDebugIntegratorStages &&
            t >= debugStartTime && t <= debugStopTime);
}
