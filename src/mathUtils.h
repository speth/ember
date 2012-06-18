#pragma once

#define USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <map>

#include <Eigen/Dense>

using std::abs;

typedef Eigen::ArrayXXd dmatrix;
typedef Eigen::ArrayXd dvec;

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

    int nanloc(const dvector& v); // returns index of first NaN component. Returns -1 if none

    dvector abs(const dvector& v);

    // Returns the index of the first/last element of v which is true
    // Returns -1 if all elements of v are false
    int findFirst(const vector<bool>& v);
    int findLast(const vector<bool>& v);

    // Returns the indices of the elements of v which are true
    vector<int> find(vector<bool>& v);

    void smooth(dvector& v);
    void smooth(dvec& v);
    dvector linspace(const double x1, const double x2, const int n);

    // Internal function for splines
    vector<dvector> computeSplines(const dvector& xIn, const dvector& yIn);

    // Linear interpolation
    dvector interp1(const dvector& xIn, const dvector& yIn, const dvector& xOut);
    double interp1(const dvector& xIn, const dvector& yIn, const double xOut);

    // Cubic Spline interpolation
    dvector splines(const dvector& xIn, const dvector& yIn, const dvector& xOut);
    double splines(const dvector& xIn, const dvector& yIn, const double xOut);

    // Integration of cubic splines
    double integrate(const dvector& x, const dvector& y);

    // Trapezoidal rule integration
    double trapz(const dvector& x, const dvector& y);

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

