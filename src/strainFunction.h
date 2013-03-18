#pragma once

class ConfigOptions;

class StrainFunction
{
public:
    StrainFunction() {};
    virtual ~StrainFunction() {}
    virtual double a(double t) const = 0;
    virtual double dadt(double t) const = 0;
};

class LinearStrainFunction : public StrainFunction
{
public:
    // Strain rate is constant at aInitial until time T0,
    // then increases linearly to aFinal at time T0+Dt
    explicit LinearStrainFunction(const ConfigOptions& options);
    virtual double a(double t) const;
    virtual double dadt(double t) const;

private:
    double aInitial, aFinal;
    double T0, Dt;
};
