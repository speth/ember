#pragma once

#define _USE_MATH_DEFINES
#define NOMINMAX
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <map>

#include "debugUtils.h"
#include "cantera/base/ct_defs.h"

#include <eigen3/Eigen/Dense>

using std::abs;
using std::size_t;
using Cantera::npos;

typedef Eigen::ArrayXXd dmatrix;
typedef Eigen::ArrayXd dvec;
typedef Eigen::ArrayXd::Index index_t;
typedef Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> StrideXX;
typedef Eigen::Stride<1,Eigen::Dynamic> Stride1X;
typedef Eigen::Map<dmatrix, Eigen::Aligned, StrideXX > MatrixMap;
typedef Eigen::Map<dvec, Eigen::Unaligned, Stride1X> VecMap;

using std::vector;
typedef std::vector<double> dvector;

class sdVector;

inline void remap(dmatrix& M, MatrixMap& A,
                  index_t nRows, index_t nCols, index_t start)
{
    assert(start + nRows <= M.rows());
    assert(nCols <= M.cols());
    new (&A) MatrixMap(&M.row(start)[0], nRows, nCols,
                       StrideXX(M.outerStride(), M.innerStride()));
}

inline void remap(dmatrix& M, VecMap& A,
                  index_t nVals, index_t start)
{
    assert(start < M.rows());
    assert(nVals <= M.cols());
    new (&A) VecMap(&M.row(start)[0], nVals, Stride1X(1, M.outerStride()));
}

extern const double NaN;

//! A collection of general-purpose mathematical functions for working with
//! vectors.
namespace mathUtils
{
    //! Returns the value of the maximum element of a vector
    double maxval(const dvector& v);

    //! Returns the value of the minimum element of a vector
    double minval(const dvector& v);

    //! Returns the range (maximum - minimum) of a vector
    double range(const dvector& v);

    //! Returns the sum of all the elements in a vector
    double sum(const dvector& v);

    //! Returns the arithemetic mean of the elements in a vector
    double mean(const dvector& v);

    //! Returns the value of the maximum element of the sub-vector of *v* with
    //! indices on the closed interval `[iStart, iEnd]`
    double maxval(const dvector& v, size_t iStart, size_t iEnd);

    //! Returns the value of the minimum element of the sub-vector of *v* with
    //! indices on the closed interval `[iStart, iEnd]`
    double minval(const dvector& v, size_t iStart, size_t iEnd);

    //! Returns the range (maximum - minimum) of the sub-vector of *v* with
    //! indices on the closed interval `[iStart, iEnd]`
    double range(const dvector& v, size_t iStart, size_t iEnd);

    //! Returns the sum of all the elements in the sub-vector of *v* with
    //! indices on the closed interval `[iStart, iEnd]`
    double sum(const dvector& v, size_t iStart, size_t iEnd);

    //! Returns the arithemetic mean of the sub-vector of *v* with indices on
    //! the closed interval `[iStart, iEnd]`
    double mean(const dvector& v, size_t iStart, size_t iEnd);

    //! Returns the index of the minimum element of *v*.
    size_t minloc(const dvector& v);

    //! Returns the index of the maximum element of *v*.
    size_t maxloc(const dvector& v);

    //! Returns true if `v` does not contain any `NaN`s
    template <typename Derived>
    bool notnan(const Eigen::EigenBase<Derived>& v)
    {
        double s = v.derived().sum();
        return (s > 0 || s <= 0);
    }

    //! Returns true if `v` does not contain any `NaN`s
    template <typename T>
    bool notnan(const vector<T>& v) {
        for (size_t i=0; i<v.size(); i++) {
            if (!(v[i] > 0) && !(v[i] <= 0)) {
                return false;
            }
        }
        return true;
    }

    //! Returns true if `v` is not `NaN`.
    inline bool notnan(const double& v) {
        return (v > 0 || v <= 0);
    }

    //! Returns true if `v` does not contain any `NaN`s
    bool notnan(const sdVector& v);

    //! Returns the index of first `NaN` component of `v`. Returns -1 if no
    //! elements are `NaN`.
    size_t nanloc(const dvector& v);

    //! Returns `true` if *a* and *b* are equal to within the specified
    //! relative (`rtol`) and absolute (`atol`) tolerances.
    inline bool almostEqual(double a, double b, double rtol=1e-10, double atol=1e-18) {
        return (abs(a - b) < rtol * (abs(a) + abs(b)) + atol);
    }

    //! Returns a vector containing the absolute values of the elements of
    //! *v*.
    dvector abs(const dvector& v);

    //! Returns the index of the first element of *v* which is `true`. Returns
    //! `-1` if all elements of *v* are `false`.
    template <class T>
    int findFirst(const T& v) {
        for (long int i=0; i< static_cast<long int>(v.size()); i++) {
            if (v[i]) {
                return i;
            }
        }
        return -1;
    }

    //! Returns the index of the last element of *v* which is `true`. Returns `-1` if
    //! all elements of *v* are `false`.
    template <class T>
    int findLast(const T& v)
    {
        for (long int i = static_cast<long int>(v.size())-1; (i+1)>0; i--) {
            if (v[i]) {
                return i;
            }
        }
        return -1;
    }

    //! Returns a vector containing the indices of the elements of *v* which
    //! are `true`.
    vector<int> find(vector<bool>& v);

    //! Applies a simple smoothing function to *v* which effectively filters
    //! out high-frequency components.
    template <class T>
    void smooth(T& v)
    {
        if (v.size() == 0) {return;}
        double p = v[0];
        double q = 0;
        for (long int i=1; i<static_cast<long int>(v.size())-1; i++) {
            q = v[i];
            v[i] = 0.25*p + 0.5*v[i] + 0.25*v[i+1];
            std::swap(p,q);
        }
    }

    //! Returns a vector of *n* elements linearly spaced along the closed
    //! interval `[x1, x2]`
    dvector linspace(const double x1, const double x2, const int n);

    //! Internal function for calculating cubic splines for a sequence of data
    //! points.
    template <class T1, class T2>
    vector<T1> computeSplines(const T1& xIn, const T2& yIn)
    {
        if (xIn.size() != yIn.size()) {
            throw DebugException("mathUtils::ComputeSplines: error: xIn and yIn must be the same size.");
        }

        int nIn = static_cast<int>(xIn.size());

        // Compute spacing of x and derivative of y
        T1 h(nIn-1), b(nIn-1);
        for (int i=0; i<nIn-1;  i++) {
            h[i] = xIn[i+1]-xIn[i];
            b[i] = (yIn[i+1]-yIn[i])/h[i];
        }

        // Gaussian elimination for the Tridiagonal system
        T1 u(nIn-1), v(nIn-1);
        u[0] = 0; v[0] = 0;
        u[1] = 2*(h[0]+h[1]);
        v[1] = 6*(b[1]-b[0]);
        for (int i=2; i<nIn-1; i++) {
            u[i] = 2*(h[i-1]+h[i]) - h[i-1]*h[i-1]/u[i-1];
            v[i] = 6*(b[i]-b[i-1]) - h[i-1]*v[i-1]/u[i-1];
        }

        // Back-substitute
        T1 z(nIn);
        z[nIn-1] = 0;
        for (int i=nIn-2; i>0; i--) {
            z[i] = (v[i]-h[i]*z[i+1])/u[i];
        }
        z[0] = 0;
        assert(notnan(z));

        // Evaluate the polynomial coefficients
        vector<T1> c(4, T1(nIn-1));
        for (int i=0; i<nIn-1; i++) {
            c[0][i] = yIn[i];
            c[1][i] = -h[i]/6*z[i+1] - h[i]/3*z[i] + (yIn[i+1]-yIn[i])/h[i];
            c[2][i] = 0.5*z[i];
            c[3][i] = (z[i+1]-z[i])/(6*h[i]);
        }
        assert(notnan(c[0]));
        assert(notnan(c[1]));
        assert(notnan(c[2]));
        assert(notnan(c[3]));

        return c;
    }

    //! Linear interpolation for a vector of desired outputs. `T` is a
    //! container of doubles, e.g. `vector<double>` or `Eigen::ArrayXd`.
    //! @param xIn Input x-coordinates, length *N*. Must be monotonically
    //!     increasing.
    //! @param yIn Input y-coordinates, length *N*
    //! @param xOut Target x-coordinates, length *M*
    //! @param extrap Linearly extrapolate outside the input interval. If false,
    //!     clip to the nearest value.
    //! @returns Computed y-coordinates at the points defined in *xOut*,
    //!     length *M*
    template <class T>
    T interp1(const T& xIn, const T& yIn, const T& xOut, bool extrap=true)
    {
        if (xIn.size() != yIn.size()) {
            throw DebugException("mathUtils::interp1: error: xIn and yIn must be the same size.");
        }
        int nOut = static_cast<int>(xOut.size());
        int nIn = static_cast<int>(xIn.size());

        T yOut(nOut);

        if (extrap) {
            for (int i=0; i<nOut; i++) {
                int j = findFirst(xIn >= xOut[i]) - 1;
                if (j == -1) {
                    j = 0;
                } else if (j == -2) {
                    j = nIn-2;
                }
                yOut[i] = yIn[j] + (yIn[j+1]-yIn[j])/(xIn[j+1]-xIn[j])*(xOut[i]-xIn[j]);
            }
        } else {
            for (int i=0; i<nOut; i++) {
                int j = findFirst(xIn >= xOut[i]) - 1;
                if (j == -1) {
                    yOut[i] = yIn[0];
                } else if (j == -2) {
                    yOut[i] = yIn[nIn-1];
                } else {
                    yOut[i] = yIn[j] + (yIn[j+1]-yIn[j])/(xIn[j+1]-xIn[j])*(xOut[i]-xIn[j]);
                }
            }
        }

        return yOut;
    }

    //! Linear interpolation for a single output point. `T` is a container of
    //! dobules, e.g. `vector<double>` or `Eigen::ArrayXd`.
    //! @param xIn Input x-coordinates, length *N*. Must be monotonically
    //!     increasing.
    //! @param yIn Input y-coordinates, length *N*
    //! @param xOut Target x-coordinate
    //! @param extrap Linearly extrapolate outside the input interval. If false,
    //!     clip to the nearest value.
    //! @returns Computed y-coordinate at *xOut*
    template <class T>
    double interp1(const T& xIn, const T& yIn, const double xOut,
                   bool extrap=true)
    {
        if (xIn.size() != yIn.size()) {
            throw DebugException("mathUtils::interp1: error: xIn and yIn must be the same size.");
        }

        int nIn = static_cast<int>(xIn.size());

        int j = findFirst(xIn >= xOut) - 1;
        if (j == -1) {
            if (!extrap) {
                return yIn[0];
            }
            j = 0;
        } else if (j == -2) {
            if (!extrap) {
                return yIn[nIn-1];
            }
            j = nIn-2;
        }
        return yIn[j] + (yIn[j+1]-yIn[j])/(xIn[j+1]-xIn[j])*(xOut-xIn[j]);
    }

    //! Cubic spline interpolation for a vector of desired outputs.
    //! @param xIn Input x-coordinates, length *N*. Must be monotonically
    //!     increasing.
    //! @param yIn Input y-coordinates, length *N*
    //! @param xOut Target x-coordinates, length *M*
    //! @returns Computed y-coordinates at the points defined in *xOut*,
    //!     length *M*
    dvec splines(const dvec& xIn, const dvec& yIn, const dvec& xOut);

    //! Cubic spline interpolation for a single output point.
    //! @param xIn Input x-coordinates, length *N*. Must be monotonically
    //!     increasing.
    //! @param yIn Input y-coordinates, length *N*
    //! @param xOut Target x-coordinate
    //! @returns Computed y-coordinate at *xOut*
    double splines(const dvec& xIn, const dvec& yIn, const double xOut);

    //! Numerical integration using the cubic spline fit for the given data.
    template <class T1, class T2>
    double integrate(const T1& x, const T2& y)
    {
        if (x.size() != y.size()) {
            throw DebugException("mathUtils::integrate: error: xIn and yIn must be the same size.");
        }

        assert(notnan(x));
        assert(notnan(y));
        int n = static_cast<int>(x.size());
        vector<T1> c = computeSplines(x, y);
        double I = 0;

        // Integrate the spline on each interval:

        for (int i=0; i<n-1; i++) {
            double dx = x[i+1]-x[i];
            I += dx*(c[0][i] + dx*(c[1][i]/2 + dx*(c[2][i]/3 + dx*c[3][i]/4)));
        }

        return I;
    }

    //! Trapezoidal rule integration
    template <class T1, class T2>
    double trapz(const T1& x, const T2& y)
    {
        if (x.size() != y.size()) {
            throw DebugException("mathUtils::trapz: error: x and y must be the same size.");
        }

        double I = 0;
        int n = static_cast<int>(x.size());
        for (int i=0; i<n-1; i++) {
            I += 0.5*(y[i+1]+y[i])*(x[i+1]-x[i]);
        }
        return I;
    }

    void vectorVectorToArray2D(const vector<dvector>& v, dmatrix& a);
    void array2DToVectorVector(const dmatrix& a, vector<dvector>& v);

    // For converting numbers into strings
    std::string stringify(double x);
    std::string stringify(double x, int nDigits);
    std::string stringify(int x);

    // Sign: why this isn't defined in cmath, I'll never know
    int sign(const double x);
    int sign(const int x);

    //! Sort & remove duplicate entries from *keys*, and perform the same
    //! permutation on *values* Requires that `keys.size() == values[i].size()`
    template <class Tx, class Ty>
    void uniqueSort(vector<Tx>& keys, vector< vector<Ty> >& values);
}

std::ostream& operator<<(std::ostream& os, const dvector& v);
std::ostream& operator<<(std::ostream& os, const vector<bool>& v);
std::ostream& operator<<(std::ostream& os, const vector<int>& v);

dvector& operator+=(dvector& v1, const dvector& v2);
dvector operator+(const dvector& v1, const dvector& v2);
dvector& operator-=(dvector& v1, const dvector& v2);
dvector operator-(const dvector& v1, const dvector& v2);

dvector& operator+=(dvector& v1, const double s);
dvector operator+(const dvector& v1, const double s);
dvector& operator-=(dvector& v1, const double s);
dvector operator-(const dvector& v1, const double s);

dvector& operator*=(dvector& v1, const dvector& v2);
dvector operator*(const dvector& v1, const dvector& v2);
dvector& operator*=(dvector& v1, const double s);
dvector operator*(const dvector& v1, const double s);
dvector operator*(const double s, const dvector& v1);

dvector& operator/=(dvector& v1, const dvector& v2);
dvector operator/(const dvector& v1, const dvector& v2);
dvector& operator/=(dvector& v1, const double s);
dvector operator/(const dvector& v1, const double s);

vector<bool> operator>(const dvector& v1, const dvector& v2);
vector<bool> operator<(const dvector& v1, const dvector& v2);
vector<bool> operator>=(const dvector& v1, const dvector& v2);
vector<bool> operator<=(const dvector& v1, const dvector& v2);
vector<bool> operator==(const dvector& v1, const dvector& v2);
vector<bool> operator!=(const dvector& v1, const dvector& v2);

vector<bool> operator>(const dvector& v1, const double& s);
vector<bool> operator<(const dvector& v1, const double& s);
vector<bool> operator>=(const dvector& v1, const double& s);
vector<bool> operator<=(const dvector& v1, const double& s);
vector<bool> operator==(const dvector& v1, const double& s);
vector<bool> operator!=(const dvector& v1, const double& s);

vector<bool> operator!(const vector<bool>& v);
vector<bool> operator&&(const vector<bool>& v1, const vector<bool>& v2);
vector<bool> operator||(const vector<bool>& v1, const vector<bool>& v2);


template <class Tx, class Ty>
void mathUtils::uniqueSort(vector<Tx>& keys, vector< vector<Ty> >& values)
{
    int n = keys.size();
    int m = values.size();
    std::map<Tx, vector<Ty> > box;
    vector<double> tmp(m);

    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            tmp[j] = values[j][i];
        }
        box[keys[i]] = tmp;
    }

    n = box.size();
    keys.resize(n);
    for (int j=0; j<m; j++) {
        values[j].resize(n);
    }

    typename std::map<Tx, vector<Ty> >::iterator iter;
    int i=0;
    for (iter=box.begin(); iter!=box.end(); iter++) {
        keys[i] = iter->first;
        for (int j=0; j<m; j++) {
            values[j][i] = iter->second[j];
        }
        i++;
    }
}

