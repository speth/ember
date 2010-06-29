#include "strainFunction.h"

StrainFunction::StrainFunction(const configOptions& options)
    : aInitial(options.strainRateInitial)
    , aFinal(options.strainRateFinal)
    , T0(options.strainRateT0)
    , Dt(options.strainRateDt)
{
    aPrev = a(T0);
    tPrev = T0;
}

double StrainFunction::a(double t) const
{
    return (t <= T0) ? aInitial
        :  (t >= T0+Dt) ? aFinal
        : aInitial + (aFinal-aInitial)*(t-T0)/Dt;
}

double StrainFunction::dadt(double t) const
{
    return (t>tPrev) ? (a(t)-aPrev)/(t-tPrev) : 0;
}

void StrainFunction::pin(double t)
{
    aPrev = a(t);
    tPrev = t;
}
