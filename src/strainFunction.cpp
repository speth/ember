#include "strainFunction.h"

double StrainFunction::a(double t) const
{
    // TODO: implement configurable strain rate
    return 100;
}

double StrainFunction::dadt(double t) const
{
    // TODO: implement configurable strain rate
    return 0;
}
