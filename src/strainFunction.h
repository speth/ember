#pragma once

#include "readConfig.h"

class StrainFunction
{
    // Strain rate is constant at aInitial until time T0,
    // then increases linearly to aFinal at time T0+Dt
public:
    StrainFunction();
    StrainFunction(const configOptions& options);

    void setOptions(const configOptions& options);
    double a(double t) const;
    double dadt(double t) const;
    void pin(double t); // call this after each successful timestep

private:
    double aInitial, aFinal;
    double T0, Dt;
    double aPrev, tPrev; // reference values for calculating dadt
};
