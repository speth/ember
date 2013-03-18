#include "strainFunction.h"
#include "readConfig.h"

LinearStrainFunction::LinearStrainFunction(const ConfigOptions& options)
    : aInitial(options.strainRateInitial)
    , aFinal(options.strainRateFinal)
    , T0(options.strainRateT0)
    , Dt(options.strainRateDt)
{
}

double LinearStrainFunction::a(double t) const
{
    return (t <= T0) ? aInitial
        :  (t >= T0+Dt) ? aFinal
        : aInitial + (aFinal-aInitial)*(t-T0)/Dt;
}

double LinearStrainFunction::dadt(double t) const
{
    return (t >= T0 && t < T0 + Dt) ? (aFinal-aInitial) / Dt : 0;
}
