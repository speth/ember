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

ChebyshevStrainFunction::ChebyshevStrainFunction(const ConfigOptions& options)
    : T0(options.strainRateT0)
    , T1(options.strainRateDt)
{
    coeffs.resize(2);
    coeffs << options.strainRateInitial, 0.0;
}

void ChebyshevStrainFunction::setCoefficients(int n, const double* data)
{
    assert(n >= 2);
    T0 = data[0];
    T1 = data[1];

    // To avoid special cases when evaluating a(t) and dadt(t), fill in zero
    // coefficients so that the there are always at least two terms in the
    // approximation.
    if (n >= 4) {
        coeffs = VecMap(const_cast<double*>(&data[2]), n-2, Stride1X(1, 1));
    } else if (n == 2) {
        coeffs.setConstant(2, 0.0);
    } else if (n == 3) {
        coeffs.setConstant(2, 0.0);
        coeffs[0] = data[2];
    }
}

double ChebyshevStrainFunction::a(double t) const
{
    // Map t onto the interval [-1,1] used to define the Chebyshev polynomial
    double x = 2 * (t - T0) / (T1 - T0) - 1;
    assert(x > -1.001 && x < 1.001);

    // Use the recurrence relation for Chebyshev polynomials of the first kind
    // to evaluate each term in the sum.
    double Cnm1 = 1; // T_0(x)
    double Cn = x; // T_1(x)
    double y = coeffs[0] + x * coeffs[1];
    for (index_t i = 2; i < coeffs.size(); i++) {
        double Cnp1 = 2 * x * Cn - Cnm1; // T_i(x)
        y += Cnp1 * coeffs[i];
        Cnm1 = Cn;
        Cn = Cnp1;
    }
    return y;
}

double ChebyshevStrainFunction::dadt(double t) const
{
    double x = 2 * (t - T0) / (T1 - T0) - 1;
    assert(x > -1.001 && x < 1.001);

    // Evaluate using the fact that T'_n(x) = n*U_(n-1)
    double Cnm1 = 1; // U_0(x)
    double Cn = 2 * x; // U_1(x)
    double dydx = coeffs[1];
    for (index_t i = 2; i < coeffs.size(); i++) {
        double Cnp1 = 2 * x * Cn - Cnm1; // U_i(x)
        dydx += i * Cnp1 * coeffs[i];
        Cnm1 = Cn;
        Cn = Cnp1;
    }

    double dxdt = 2.0 / (T1 - T0);
    return dydx * dxdt;
}

StrainFunction* newStrainFunction(const ConfigOptions& options)
{
    if (options.strainFunctionType == "linear") {
        return new LinearStrainFunction(options);
    } else if (options.strainFunctionType == "chebyshev") {
        return new ChebyshevStrainFunction(options);
    } else {
        throw DebugException("No such StrainFunction type: '" +
                             options.strainFunctionType + "'.");
    }
}
