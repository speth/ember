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

#include <Eigen/Dense>

using std::abs;
using std::size_t;

typedef Eigen::ArrayXXd dmatrix;
typedef Eigen::ArrayXd dvec;
typedef Eigen::ArrayXd::Index index_t;

using std::vector;
typedef std::vector<double> dvector;

extern const double NaN;

namespace mathUtils
{
    double maxval(const dvector& v);
    double minval(const dvector& v);
    double range(const dvector& v);
    double sum(const dvector& v);
    double mean(const dvector& v);

    // Find the min/max of a subvector
    double maxval(const dvector& v, int iStart, int iEnd);
    double minval(const dvector& v, int iStart, int iEnd);
    double range(const dvector& v, int iStart, int iEnd);
    double sum(const dvector& v, int iStart, int iEnd);
    double mean(const dvector& v, int iStart, int iEnd);

    // location of the minimum / maximum
    int minloc(const dvector& v);
    int maxloc(const dvector& v);

    //! Returns true if v does not contain any NaNs
    template <class T>
    bool notnan(const T& v) {
        for (size_t i=0; i<v.size(); i++) {
            if (!(v[i] > 0) && !(v[i] <= 0)) {
                return false;
            }
        }
        return true;
    }

    template <>
    bool notnan<dmatrix>(const dmatrix& v);

    template <>
    inline bool notnan<double>(const double& v) {
        return (v > 0 || v <= 0);
    }

    int nanloc(const dvector& v); // returns index of first NaN component. Returns -1 if none

    dvector abs(const dvector& v);

    // Returns the index of the first element of v which is true
    // Returns -1 if all elements of v are false.
    template <class T>
    int findFirst(const T& v) {
        for (long int i=0; i< static_cast<long int>(v.size()); i++) {
            if (v[i]) {
                return i;
            }
        }
        return -1;
    }

    // Returns the index of the last element of v which is true
    // Returns -1 if all elements of v are false.
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

    // Returns the indices of the elements of v which are true
    vector<int> find(vector<bool>& v);

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

    dvector linspace(const double x1, const double x2, const int n);

    // Internal function for splines
    template <class T1, class T2>
    vector<T1> computeSplines(const T1& xIn, const T2& yIn)
    {
        if (xIn.size() != yIn.size()) {
            throw debugException("mathUtils::ComputeSplines: error: xIn and yIn must be the same size.");
        }

        int nIn = xIn.size();

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

    // Linear interpolation. T is, e.g. vector<double>.
    template <class T>
    T interp1(const T& xIn, const T& yIn, const T& xOut)
    {
        if (xIn.size() != yIn.size()) {
            throw debugException("mathUtils::interp1: error: xIn and yIn must be the same size.");
        }
        int nOut = xOut.size();
        int nIn = xIn.size();

        T yOut(nOut);

        for (int i=0; i<nOut; i++) {
            int j = findFirst(xIn >= xOut[i]) - 1;
            if (j == -1) {
                j = 0;
            } else if (j == -2) {
                j = nIn-2;
            }
            yOut[i] = yIn[j] + (yIn[j+1]-yIn[j])/(xIn[j+1]-xIn[j])*(xOut[i]-xIn[j]);
        }

        return yOut;
    }

    // Linear interpolation. T is a container of doubles.
    template <class T>
    double interp1(const T& xIn, const T& yIn, const double xOut)
    {
        if (xIn.size() != yIn.size()) {
            throw debugException("mathUtils::interp1: error: xIn and yIn must be the same size.");
        }

        int nIn = xIn.size();

        int j = findFirst(xIn >= xOut) - 1;
        if (j == -1) {
            j = 0;
        } else if (j == -2) {
            j = nIn-2;
        }
        return yIn[j] + (yIn[j+1]-yIn[j])/(xIn[j+1]-xIn[j])*(xOut-xIn[j]);
    }

    // Cubic Spline interpolation
    dvec splines(const dvec& xIn, const dvec& yIn, const dvec& xOut);
    double splines(const dvec& xIn, const dvec& yIn, const double xOut);

    // Integration of cubic splines
    template <class T1, class T2>
    double integrate(const T1& x, const T2& y)
    {
        if (x.size() != y.size()) {
            throw debugException("mathUtils::integrate: error: xIn and yIn must be the same size.");
        }

        assert(notnan(x));
        assert(notnan(y));
        int n = x.size();
        vector<T1> c = computeSplines(x, y);
        double I = 0;

        // Integrate the spline on each interval:

        for (int i=0; i<n-1; i++) {
            double dx = x[i+1]-x[i];
            I += dx*(c[0][i] + dx*(c[1][i]/2 + dx*(c[2][i]/3 + dx*c[3][i]/4)));
        }

        return I;
    }

    // Trapezoidal rule integration
    template <class T1, class T2>
    double trapz(const T1& x, const T2& y)
    {
        if (x.size() != y.size()) {
            throw debugException("mathUtils::trapz: error: x and y must be the same size.");
        }

        double I = 0;
        int n = x.size();
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

    // Sort & remove duplicate entries from "keys", and perform the same permutation on "values"
    // Requires that keys.size() == values[i].size()
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

