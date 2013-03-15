#pragma once

#include "mathUtils.h"

class ConfigOptions;

class StrainFunction
{
public:
    StrainFunction() {};
    virtual ~StrainFunction() {}
    virtual void setCoefficients(int n, const double* coeffs) {}

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

//! Strain rate calculated as a Chebyshev polynomial representation of a(t).
/*! \f[
 *      a(t) = \sum_{n=0}^N a_n T_n(x)
 *  \f]
 *  where \f$ T_n(x) \f$ are Chebyshev polynomials of the first kind and
 * \f$ x = 2 * (t-t_o)/(t_1-t_0) - 1 \f$ maps *t* onto the interval [-1,1].
 */
class ChebyshevStrainFunction: public StrainFunction
{
public:
    explicit ChebyshevStrainFunction(const ConfigOptions& options);

    virtual double a(double t) const;
    virtual double dadt(double t) const;

    //! The first two elements are the time interval \f$ (t_0, t_1) \f$ used
    //! to map the time onto the interval (-1, 1). The remaining elements are
    //! the coefficients used to multiply the Chebyshev polynomials.
    virtual void setCoefficients(int n, const double* coeffs);

private:
    double T0, T1;
    dvec coeffs;
};

StrainFunction* newStrainFunction(const ConfigOptions& options);
