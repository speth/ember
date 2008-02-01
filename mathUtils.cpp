#include "mathUtils.h"

double mathUtils::max(const vector<double>& v)
{
	//vector<double>::const_iterator i;
	//double val = *v.begin();
	//for (i=v.begin(); i!=v.end(); i++) {
	//	val = std::max(*i,val);
	//}
	//return val;
	int n = v.size();
	if (n==0) { return 0; }
	double val = v[0];
	for (int i=0; i<n; i++) {
		if (v[i] > val) {
			val = v[i];
		}
	}
	return val;
}

double mathUtils::min(vector<double>& v)
{
	vector<double>::const_iterator i;
	double val = *v.begin();
	for (i=v.begin(); i!=v.end(); i++) {
		val = std::min(*i,val);
	}
	return val;
}

double mathUtils::range(vector<double>& v)
{
	return (max(v)-min(v));
}

int mathUtils::minloc(vector<bool>& v)
{
	for (vector<bool>::size_type i=0; i<v.size(); i++) {
		if (v[i]) {
			return i;
		}
	}
	return -1;
}


int mathUtils::maxloc(vector<bool>& v)
{
	for (vector<bool>::size_type i=v.size()-1; i>=0; i--) {
		if (v[i]) {
			return i;
		}
	}
	return -1;
}

vector<int> mathUtils::find(vector<bool>& v)
{
	vector<int> out;
	for (vector<bool>::size_type i=0; i<v.size(); i++) {
		if (v[i]) {
			out.push_back(i);
		}
	}
	return out;
}

std::ostream& operator<<(std::ostream& os, vector<double>& v)
{
	for (vector<double>::size_type i=0; i<v.size()-1; i++) {
		os << v[i] << ", ";
	}
	os << v[v.size()-1];
	return os;
}

std::ostream& operator<<(std::ostream& os, vector<bool>& v)
{
	for (vector<bool>::size_type i=0; i<v.size()-1; i++) {
		os << v[i] << ", ";
	}
	os << v[v.size()-1];
	return os;
}

std::ostream& operator<<(std::ostream& os, vector<int>& v)
{
	for (vector<int>::size_type i=0; i<v.size()-1; i++) {
		os << v[i] << ", ";
	}
	os << v[v.size()-1];
	return os;
}

vector<double>& operator+=(vector<double>& v1, const vector<double>& v2)
{
	if (v1.size() != v2.size()) {
		throw;
	}
	unsigned int n = v1.size();

	for (unsigned int i=0; i<n; i++) {
		v1[i] += v2[i];
	}

	return v1;
}

vector<double> operator+(const vector<double>& v1, const vector<double>& v2)
{
	vector<double> v(v1);
	return (v += v2);
}

vector<bool> operator>(const dvector& v1, const dvector& v2)
{
	if (v1.size() != v2.size()) throw;
	dvector::size_type n = v1.size();
	vector<bool> v(n);
	for (dvector::size_type i=0; i<n; i++) {
		v[i] = v1[i] > v2[i];
	}
	return v;
}

vector<bool> operator<(const dvector& v1, const dvector& v2)
{
	if (v1.size() != v2.size()) throw;
	dvector::size_type n = v1.size();
	vector<bool> v(n);
	for (dvector::size_type i=0; i<n; i++) {
		v[i] = v1[i] < v2[i];
	}
	return v;
}

vector<bool> operator>=(const dvector& v1, const dvector& v2)
{
	if (v1.size() != v2.size()) throw;
	dvector::size_type n = v1.size();
	vector<bool> v(n);
	for (dvector::size_type i=0; i<n; i++) {
		v[i] = v1[i] >= v2[i];
	}
	return v;
}

vector<bool> operator<=(const dvector& v1, const dvector& v2)
{
	if (v1.size() != v2.size()) throw;
	dvector::size_type n = v1.size();
	vector<bool> v(n);
	for (dvector::size_type i=0; i<n; i++) {
		v[i] = v1[i] <= v2[i];
	}
	return v;
}

vector<bool> operator==(const dvector& v1, const dvector& v2)
{
	if (v1.size() != v2.size()) throw;
	dvector::size_type n = v1.size();
	vector<bool> v(n);
	for (dvector::size_type i=0; i<n; i++) {
		v[i] = v1[i] == v2[i];
	}
	return v;
}

vector<bool> operator!=(const dvector& v1, const dvector& v2)
{
	if (v1.size() != v2.size()) throw;
	dvector::size_type n = v1.size();
	vector<bool> v(n);
	for (dvector::size_type i=0; i<n; i++) {
		v[i] = v1[i] != v2[i];
	}
	return v;
}

vector<bool> operator>(const dvector& v1, const double& s)
{
	dvector::size_type n = v1.size();
	vector<bool> v(n);
	for (dvector::size_type i=0; i<n; i++) {
		v[i] = v1[i] > s;
	}
	return v;
}

vector<bool> operator<(const dvector& v1, const double& s)
{
	dvector::size_type n = v1.size();
	vector<bool> v(n);
	for (dvector::size_type i=0; i<n; i++) {
		v[i] = v1[i] < s;
	}
	return v;
}

vector<bool> operator>=(const dvector& v1, const double& s)
{
	dvector::size_type n = v1.size();
	vector<bool> v(n);
	for (dvector::size_type i=0; i<n; i++) {
		v[i] = v1[i] >= s;
	}
	return v;
}

vector<bool> operator<=(const dvector& v1, const double& s)
{
	dvector::size_type n = v1.size();
	vector<bool> v(n);
	for (dvector::size_type i=0; i<n; i++) {
		v[i] = v1[i] <= s;
	}
	return v;
}

vector<bool> operator==(const dvector& v1, const double& s)
{
	dvector::size_type n = v1.size();
	vector<bool> v(n);
	for (dvector::size_type i=0; i<n; i++) {
		v[i] = v1[i] == s;
	}
	return v;
}

vector<bool> operator!=(const dvector& v1, const double& s)
{
	dvector::size_type n = v1.size();
	vector<bool> v(n);
	for (dvector::size_type i=0; i<n; i++) {
		v[i] = v1[i] != s;
	}
	return v;
}

vector<bool> operator!(const vector<bool>& v)
{
	vector<bool>::size_type n = v.size();
	vector<bool> out(n);
	for (vector<bool>::size_type i=0; i<n; i++) {
		out[i] = !v[i];
	}
	return out;
}

vector<bool> operator&&(const vector<bool>& v1, const vector<bool>& v2)
{
	vector<bool>::size_type n = v1.size();
	if (n != v2.size()) throw;
	vector<bool> out(n);
	for (vector<bool>::size_type i=0; i<n; i++) {
		out[i] = v1[i] && v2[i];
	}
	return out;
}

vector<bool> operator||(const vector<bool>& v1, const vector<bool>& v2)
{
	vector<bool>::size_type n = v1.size();
	if (n != v2.size()) throw;
	vector<bool> out(n);
	for (vector<bool>::size_type i=0; i<n; i++) {
		out[i] = v1[i] || v2[i];
	}
	return out;
}
