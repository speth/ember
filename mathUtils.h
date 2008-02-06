#pragma once

#include <vector>
#include <cmath>
#include <valarray>
#include <iostream>
#include <string>

using std::vector;
using std::valarray;

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

	// Returns the index of the first element of v which is true
	// Returns -1 if all elements of v are false
	int minloc(vector<bool>& v);
	int maxloc(vector<bool>& v);

	// Returns the indices of the elements of v which are true
	vector<int> find(vector<bool>& v);
};

std::ostream& operator<<(std::ostream& os, vector<double>& v);
std::ostream& operator<<(std::ostream& os, vector<bool>& v);
std::ostream& operator<<(std::ostream& os, vector<int>& v);

dvector& operator+=(dvector& v1, const dvector& v2);
dvector operator+(const dvector& v1, const dvector& v2);

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
