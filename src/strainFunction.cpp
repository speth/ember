#include "strainFunction.h"
#include "readConfig.h"

double StrainFunction::dadt(double t) const
{
    return (t>tPrev) ? (a(t)-aPrev)/(t-tPrev) : 0;
}

void StrainFunction::pin(double t)
{
    aPrev = a(t);
    tPrev = t;
}

LinearStrainFunction::LinearStrainFunction(const ConfigOptions& options)
    : aInitial(options.strainRateInitial)
    , aFinal(options.strainRateFinal)
    , T0(options.strainRateT0)
    , Dt(options.strainRateDt)
{
    pin(T0);
}

double LinearStrainFunction::a(double t) const
{
    return (t <= T0) ? aInitial
        :  (t >= T0+Dt) ? aFinal
        : aInitial + (aFinal-aInitial)*(t-T0)/Dt;
}
