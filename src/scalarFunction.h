#pragma once

#include "mathUtils.h"

class ConfigOptions;

//! Represents a parameterization of a scalar function of one variable, a(t).
class ScalarFunction
{
public:
    ScalarFunction() {};
    virtual ~ScalarFunction() {}

    //! Set the parameters associated with this function. The meaning of the
    //! elements of `coeffs` depends on the parameterization.
    virtual void setCoefficients(int n, const double* coeffs) {}

    //! Compute the value of the function at time `t`.
    virtual double a(double t) const = 0;

    //! Compute the time derivative of the function at time `t`.
    virtual double dadt(double t) const = 0;
};

//! Represents a function that is constant before and after a linear ramp.
/*!
 *  The value is constant at #aInitial until time #T0, then increases
 *  linearly to #aFinal at time #T0 + #Dt.
 */
class LinearFunction : public ScalarFunction
{
public:
    explicit LinearFunction(const ConfigOptions& options);
    virtual double a(double t) const;
    virtual double dadt(double t) const;

private:
    double aInitial, aFinal;
    double T0, Dt;
};

//! Chebyshev polynomial representation of a(t).
/*! \f[
 *      a(t) = \sum_{n=0}^N a_n T_n(x)
 *  \f]
 *  where \f$ T_n(x) \f$ are Chebyshev polynomials of the first kind and
 * \f$ x = 2 * (t-t_o)/(t_1-t_0) - 1 \f$ maps *t* onto the interval [-1,1].
 */
class ChebyshevFunction: public ScalarFunction
{
public:
    explicit ChebyshevFunction(const ConfigOptions& options);

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

//! Create a new function of `type` with additional options specified in
//! `options`.
ScalarFunction* newScalarFunction(const std::string& type,
                                  const ConfigOptions& options);
