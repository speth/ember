#pragma once

#include <vector>
#include <cmath>
#include <valarray>
#include <iostream>
#include <string>

#include <cantera/Cantera.h>
#include <cantera/kernel/Array.h>

using std::vector;
using std::valarray;
using std::abs;
using Cantera::Array2D;

typedef vector<double> dvector;

namespace mathUtils
{
	double max(const vector<double>& v);
	double min(const vector<double>& v);
	double range(const vector<double>& v);

	// Find the min/max of a subvector
	double max(const vector<double>& v, int iStart, int iEnd);
	double min(const vector<double>& v, int iStart, int iEnd);
	double range(const vector<double>& v, int iStart, int iEnd);

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

	// Cubic Spline interpolation
	dvector splines(const dvector& xIn, const dvector& yIn, const dvector& xOut);
	double splines(const dvector& xIn, const dvector& yIn, const double xOut);

	// Internal function for splines
	vector<dvector> computeSplines(const dvector& xIn, const dvector& yIn);

	void vectorVectorToArray2D(const vector<dvector>& v, Array2D& a);
	void array2DToVectorVector(const Array2D& a, vector<dvector>& v);
	
};

std::ostream& operator<<(std::ostream& os, vector<double>& v);
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
