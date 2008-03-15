#include "mathUtils.h"
#include "sundialsUtils.h"

using std::cout; using std::endl;

double mathUtils::max(const vector<double>& v)
{
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

double mathUtils::min(const vector<double>& v)
{
	vector<double>::const_iterator i;
	double val = *v.begin();
	for (i=v.begin(); i!=v.end(); i++) {
		val = std::min(*i,val);
	}
	return val;
}

double mathUtils::range(const vector<double>& v)
{
	return (max(v)-min(v));
}

double mathUtils::max(const vector<double>& v, int iStart, int iEnd)
{
	int n = v.size();
	if (n==0) { return 0; }
	if (iStart < 0 || iEnd > n) { throw; }
	double val = v[0];
	for (int i=iStart; i<iEnd; i++) {
		if (v[i] > val) {
			val = v[i];
		}
	}
	return val;
}

double mathUtils::min(const vector<double>& v, int iStart, int iEnd)
{
	int n = v.size();
	if (n==0) { return 0; }
	if (iStart < 0 || iEnd > n) { throw; }
	double val = v[0];
	for (int i=iStart; i<iEnd; i++) {
		if (v[i] < val) {
			val = v[i];
		}
	}
	return val;
}

int mathUtils::maxloc(const vector<double>& v)
{
	int n = v.size();
	if (n==0) { return 0; }
	double val = v[0];
	int loc = 0;
	for (int i=0; i<n; i++) {
		if (v[i] > val) {
			val = v[i];
			loc = i;
		}
	}
	return loc;
}

int mathUtils::minloc(const vector<double>& v)
{
	int n = v.size();
	if (n==0) { return 0; }
	double val = v[0];
	int loc = 0;
	for (int i=0; i<n; i++) {
		if (v[i] < val) {
			val = v[i];
			loc = i;
		}
	}
	return loc;
}

double mathUtils::range(const vector<double>& v, int iStart, int iEnd)
{
	return max(v,iStart,iEnd) - min(v,iStart,iEnd);
}

int mathUtils::findFirst(vector<bool>& v)
{
	for (vector<bool>::size_type i=0; i<v.size(); i++) {
		if (v[i]) {
			return i;
		}
	}
	return -1;
}

dvector mathUtils::abs(const dvector& v)
{
	dvector a(v);
	for (dvector::size_type i=0; i<v.size(); i++)
	{
		a[i] = std::abs(a[i]);
	}
	return a;
}

int mathUtils::findLast(vector<bool>& v)
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

void mathUtils::smooth(dvector& v) 
{
	if (v.size() == 0) {return;}
	double p = v[0];
	double q = 0;
	for (dvector::size_type i=1; i<v.size()-1; i++) {
		q = v[i];
		v[i] = 0.25*p + 0.5*v[i] + 0.25*v[i+1];
		std::swap(p,q);
	}
}

dvector mathUtils::linspace(const double x1, const double x2, const int n)
{
	dvector v(n);
	for (int i=0; i<n; i++) {
		v[i] = x1 + ((x2-x1)*i)/(n-1);
	}
	return v;
}

vector<dvector> mathUtils::computeSplines(const dvector& xIn, const dvector& yIn)
{
	if (xIn.size() != yIn.size()) {
		cout << "mathUtils::ComputeSplines: error: xIn and yIn must be the same size." << endl;
		throw;
	}

	int nIn = xIn.size();

	// Compute spacing of x and derivative of y
	dvector h(nIn), b(nIn);
	for (int i=0; i<nIn-1;  i++) {
		h[i] = xIn[i+1]-xIn[i];
		b[i] = (yIn[i+1]-yIn[i])/h[i];
	}

	// Gaussian elimination for the Tridiagonal system
	dvector u(nIn), v(nIn);
	u[0] = 0; v[0] = 0;
	u[1] = 2*(h[0]+h[1]);
	v[1] = 6*(b[1]-b[0]);
	for (int i=2; i<nIn-1; i++) {
		u[i] = 2*(h[i-1]+h[i]) - h[i-1]*h[i-1]/u[i-1];
		v[i] = 6*(b[i]-b[i-1]) - h[i-1]*v[i-1]/u[i-1];
	}

	// Back-substitute
	dvector z(nIn);
	z[nIn-1] = 0;
	for (int i=nIn-2; i>0; i--) {
		z[i] = (v[i]-h[i]*z[i+1])/u[i];
	}
	
	// Evaluate the polynomial coefficients
	dvector c0(nIn), c1(nIn), c2(nIn), c3(nIn);
	vector<dvector> c(4, dvector(nIn));
	for (int i=0; i<nIn-1; i++) {
		c[0][i] = yIn[i];
		c[1][i] = -h[i]/6*z[i+1] - h[i]/3*z[i] + (yIn[i+1]-yIn[i])/h[i];
		c[2][i] = 0.5*z[i];
		c[3][i] = (z[i+1]-z[i])/(6*h[i]);
	}

	return c;
}

dvector mathUtils::splines(const dvector& xIn, const dvector& yIn, const dvector& xOut)
{
	if (xIn.size() != yIn.size()) {
		cout << "mathUtils::splines: error: xIn and yIn must be the same size." << endl;
		throw;
	}

	int nOut = xOut.size();
	int nIn = xIn.size();

	vector<dvector> c = computeSplines(xIn, yIn);

	// Evaluate the spline at the selected points
	double dx;
	dvector yOut(nOut);

	for (int i=0; i<nOut; i++) {
		int j = findFirst(xIn >= xOut[i]) - 1;
		if (j == -1) {
			j = 0;
		} else if (j == -2) {
			j = nIn-1;
		}
		dx = xOut[i]-xIn[j];
		yOut[i] = c[0][j] + dx*(c[1][j] + dx*(c[2][j] + dx*c[3][j]));
	}

	return yOut;
}

double mathUtils::splines(const dvector& xIn, const dvector& yIn, const double xOut)
{
	if (xIn.size() != yIn.size()) {
		cout << "mathUtils::splines: error: xIn and yIn must be the same size." << endl;
		throw;
	}

	int nIn = xIn.size();

	vector<dvector> c = computeSplines(xIn, yIn);

	// Evaluate the spline at the selected point
	int j = findFirst(xIn >= xOut) - 1;
	if (j == -1) {
		j = 0;
	} else if (j == -2) {
		j = nIn-1;
	}
	double dx = xOut-xIn[j];
	return c[0][j] + dx*(c[1][j] + dx*(c[2][j] + dx*c[3][j]));
}

void mathUtils::vectorVectorToArray2D(const vector<dvector>& v, Array2D& a)
{
	int m = v.size();
	int n = v.begin()->size();
	a.resize(m,n);
	for (int i=0; i<m; i++) {
		for (int j=0; j<n; j++) {
			a(i,j) = v[i][j];
		}
	}
}

void mathUtils::array2DToVectorVector(const Array2D& a, vector<dvector>& v)
{
	int n = a.nRows();
	int m = a.nColumns();
	v.resize(n);
	for (int i=0; i<n; i++) {
		v[i].resize(m);
		for (int j=0; j<m; j++) {
			v[i][j] = a(i,j);
		}
	}
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

vector<double>& operator-=(vector<double>& v1, const vector<double>& v2)
{
	if (v1.size() != v2.size()) {
		throw;
	}
	unsigned int n = v1.size();

	for (unsigned int i=0; i<n; i++) {
		v1[i] -= v2[i];
	}
	return v1;
}

vector<double> operator-(const vector<double>& v1, const vector<double>& v2)
{
	vector<double> v(v1);
	return (v -= v2);
}

vector<double>& operator+=(vector<double>& v1, const double s)
{
	unsigned int n = v1.size();
	for (unsigned int i=0; i<n; i++) {
		v1[i] += s;
	}
	return v1;
}

vector<double> operator+(const vector<double>& v1, const double s)
{
	vector<double> v(v1);
	return (v += s);
}

vector<double>& operator-=(vector<double>& v1, const double s)
{
	unsigned int n = v1.size();
	for (unsigned int i=0; i<n; i++) {
		v1[i] -= s;
	}
	return v1;
}

vector<double> operator-(const vector<double>& v1, const double s)
{
	vector<double> v(v1);
	return (v -= s);
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

