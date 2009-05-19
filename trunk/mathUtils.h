#pragma once

#define USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include <cantera/Cantera.h>
#include <cantera/kernel/Array.h>

using std::abs;
using Cantera::Array2D;

using std::vector;
typedef std::vector<double> dvector;

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

	dvector abs(const dvector& v);

	// Returns the index of the first/last element of v which is true
	// Returns -1 if all elements of v are false
	int findFirst(const vector<bool>& v);
	int findLast(const vector<bool>& v);

	// Returns the indices of the elements of v which are true
	vector<int> find(vector<bool>& v);

	void smooth(dvector& v);
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

	void vectorVectorToArray2D(const vector<dvector>& v, Array2D& a);
	void array2DToVectorVector(const Array2D& a, vector<dvector>& v);

	// For converting numbers into strings
	std::string stringify(double x);
	std::string stringify(double x, int nDigits);
	std::string stringify(int x);
    
    // Sign: why this isn't defined in cmath, I'll never know
    int sign(const double x);
    int sign(const int x);
    
}

std::ostream& operator<<(std::ostream& os, dvector& v);
std::ostream& operator<<(std::ostream& os, vector<bool>& v);
std::ostream& operator<<(std::ostream& os, vector<int>& v);

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
