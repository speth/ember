#include "mathUtils.h"
#include "sundialsUtils.h"
#include "debugUtils.h"
#include <map>
#include <limits>

const double NaN = std::numeric_limits<double>::quiet_NaN();

double mathUtils::maxval(const dvector& v)
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

double mathUtils::minval(const dvector& v)
{
    double val = v[0];
    dvector::const_iterator iter;
    for (iter=v.begin(); iter!=v.end(); iter++) {
        val = std::min(*iter,val);
    }
    return val;
}

double mathUtils::range(const dvector& v)
{
    return (maxval(v)-minval(v));
}

double mathUtils::sum(const dvector& v)
{
    int n = v.size();
    if (n==0) { return 0; }
    double val = 0;
    for (int i=0; i<n; i++) {
        val += v[i];
    }
    return val;
}

double mathUtils::mean(const dvector& v)
{
    int n = v.size();
    if (n==0) { return 0; }
    double val = 0;
    for (int i=0; i<n; i++) {
        val += v[i];
    }
    return val/n;
}

double mathUtils::maxval(const dvector& v, int iStart, int iEnd)
{
    int n = v.size();
    if (n==0) { return 0; }
    if (iStart < 0 || iEnd > n)
        throw debugException("mathUtils::maxval: bad range specified");

    double val = v[iStart];
    for (int i=iStart; i<=iEnd; i++) {
        if (v[i] > val) {
            val = v[i];
        }
    }
    return val;
}

double mathUtils::minval(const dvector& v, int iStart, int iEnd)
{
    int n = v.size();
    if (n==0) { return 0; }
    if (iStart < 0 || iEnd > n)
        throw debugException("mathUtils::minval: bad range specified");

    double val = v[iStart];
    for (int i=iStart; i<=iEnd; i++) {
        if (v[i] < val) {
            val = v[i];
        }
    }
    return val;
}

double mathUtils::sum(const dvector& v, int iStart, int iEnd)
{
    int n = v.size();
    if (n==0) { return 0; }
    if (iStart < 0 || iEnd > n)
        throw debugException("mathUtils::sum: bad range specified");

    double val = 0;
    for (int i=iStart; i<=iEnd; i++) {
        val += v[i];
    }
    return val;
}

double mathUtils::mean(const dvector& v, int iStart, int iEnd)
{
    int n = v.size();
    if (n==0) { return 0; }
    if (iStart < 0 || iEnd > n)
        throw debugException("mathUtils::mean: bad range specified");

    double val = 0;
    for (int i=iStart; i<=iEnd; i++) {
        val += v[i];
    }
    return val/(iEnd-iStart+1);
}

int mathUtils::maxloc(const dvector& v)
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

int mathUtils::minloc(const dvector& v)
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

double mathUtils::range(const dvector& v, int iStart, int iEnd)
{
    return maxval(v,iStart,iEnd) - minval(v,iStart,iEnd);
}

int mathUtils::nanloc(const dvector& v) {
    for (size_t i=0; i<v.size(); i++) {
        if (!(v[i] > 0) && !(v[i] <= 0)) {
            return i;
        }
    }
    return -1;
}

int mathUtils::findFirst(const vector<bool>& v)
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

int mathUtils::findLast(const vector<bool>& v)
{
    for (vector<bool>::size_type i=v.size()-1; (i+1)>0; i--) {
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
        throw debugException("mathUtils::ComputeSplines: error: xIn and yIn must be the same size.");
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

dvector mathUtils::interp1(const dvector& xIn, const dvector& yIn, const dvector& xOut)
{
    if (xIn.size() != yIn.size()) {
        throw debugException("mathUtils::interp1: error: xIn and yIn must be the same size.");
    }
    int nOut = xOut.size();
    int nIn = xIn.size();

    dvector yOut(nOut);

    for (int i=0; i<nOut; i++) {
        int j = findFirst(xIn >= xOut[i]) - 1;
        if (j == -1) {
            j = 0;
        } else if (j == -2) {
            j = nIn-1;
        }
        yOut[i] = yIn[j] + (yIn[j+1]-yIn[j])/(xIn[j+1]-xIn[j])*(xOut[i]-xIn[j]);
    }

    return yOut;
}

double mathUtils::interp1(const dvector& xIn, const dvector& yIn, const double xOut)
{
    if (xIn.size() != yIn.size()) {
        throw debugException("mathUtils::interp1: error: xIn and yIn must be the same size.");
    }

    int nIn = xIn.size();

    int j = findFirst(xIn >= xOut) - 1;
    if (j == -1) {
        j = 0;
    } else if (j == -2) {
        j = nIn-1;
    }
    return yIn[j] + (yIn[j+1]-yIn[j])/(xIn[j+1]-xIn[j])*(xOut-xIn[j]);
}


dvector mathUtils::splines(const dvector& xIn, const dvector& yIn, const dvector& xOut)
{
    if (xIn.size() != yIn.size()) {
        throw debugException((boost::format(
            "mathUtils::splines: error: xIn (%i) and yIn (%i) must be the same size.") %
            xIn.size() % yIn.size()).str());
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
        throw debugException((boost::format(
            "mathUtils::splines: error: xIn (%i) and yIn (%i) must be the same size.") %
            xIn.size() % yIn.size()).str());
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

double mathUtils::integrate(const dvector& x, const dvector& y)
{
    if (x.size() != y.size()) {
        throw debugException("mathUtils::integrate: error: xIn and yIn must be the same size.");
    }

    int n = x.size();
    vector<dvector> c = computeSplines(x, y);
    double I=0;

    // Integrate the spline on each interval:

    for (int i=0; i<n-1; i++) {
        double dx = x[i+1]-x[i];
        I += dx*(c[0][i] + dx*(c[1][i]/2 + dx*(c[2][i]/3 + dx*c[3][i]/4)));
    }

    return I;
}

double mathUtils::trapz(const dvector& x, const dvector& y)
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

void mathUtils::vectorVectorToArray2D(const vector<dvector>& v, dmatrix& a)
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

void mathUtils::array2DToVectorVector(const dmatrix& a, vector<dvector>& v)
{
    int n = a.rows();
    int m = a.cols();
    v.resize(n);
    for (int i=0; i<n; i++) {
        v[i].resize(m);
        for (int j=0; j<m; j++) {
            v[i][j] = a(i,j);
        }
    }
}

std::string mathUtils::stringify(double x)
{
    std::ostringstream os;
    os << x;
    return os.str();
}

std::string mathUtils::stringify(int x)
{
    std::ostringstream os;
    os << x;
    return os.str();
}

std::string mathUtils::stringify(double x, int nDigits)
{
    std::ostringstream os;
    os.fill('0');
    os.width(nDigits);
    os << x;
    return os.str();
}

int mathUtils::sign(const double x)
{
    return (x > 0) ? 1 :
           (x < 0) ? -1
                   : 0;
}

int mathUtils::sign(const int x)
{
    return (x > 0) ? 1 :
           (x < 0) ? -1
                   : 0;
}

std::ostream& operator<<(std::ostream& os, const dvector& v)
{
    os << "[";
    for (dvector::size_type i=0; i<v.size()-1; i++) {
        os << v[i] << ", ";
    }
    os << v[v.size()-1];
    os << "]";
    return os;
}

std::ostream& operator<<(std::ostream& os, const vector<bool>& v)
{
    os << "[";
    for (vector<bool>::size_type i=0; i<v.size()-1; i++) {
        os << v[i] << ", ";
    }
    os << v[v.size()-1];
    os << "]";
    return os;
}

std::ostream& operator<<(std::ostream& os, const vector<int>& v)
{
    os << "[";
    for (vector<int>::size_type i=0; i<v.size()-1; i++) {
        os << v[i] << ", ";
    }
    os << v[v.size()-1];
    os << "]";
    return os;
}

dvector& operator+=(dvector& v1, const dvector& v2)
{
    if (v1.size() != v2.size())
        throw debugException("mathUtils::operator+=: cannot add vectors of different sizes.");

    unsigned int n = v1.size();

    for (unsigned int i=0; i<n; i++) {
        v1[i] += v2[i];
    }

    return v1;
}

dvector operator+(const dvector& v1, const dvector& v2)
{
    dvector v(v1);
    return (v += v2);
}

dvector& operator-=(dvector& v1, const dvector& v2)
{
    if (v1.size() != v2.size())
        throw debugException("mathUtils::operator-=: cannot subtract vectors of different sizes.");

    unsigned int n = v1.size();

    for (unsigned int i=0; i<n; i++) {
        v1[i] -= v2[i];
    }
    return v1;
}

dvector operator-(const dvector& v1, const dvector& v2)
{
    dvector v(v1);
    return (v -= v2);
}

dvector& operator+=(dvector& v1, const double s)
{
    unsigned int n = v1.size();
    for (unsigned int i=0; i<n; i++) {
        v1[i] += s;
    }
    return v1;
}

dvector operator+(const dvector& v1, const double s)
{
    dvector v(v1);
    return (v += s);
}

dvector& operator-=(dvector& v1, const double s)
{
    unsigned int n = v1.size();
    for (unsigned int i=0; i<n; i++) {
        v1[i] -= s;
    }
    return v1;
}

dvector operator-(const dvector& v1, const double s)
{
    dvector v(v1);
    return (v -= s);
}

dvector& operator*=(dvector& v1, const double s)
{
    unsigned int n = v1.size();
    for (unsigned int i=0; i<n; i++) {
        v1[i] *= s;
    }
    return v1;
}

dvector& operator*=(dvector& v1, const dvector& v2)
{
    if (v1.size() != v2.size()) {
        throw debugException("mathUtils::operator*=: cannot multiply vectors of different sizes.");
    }
    unsigned int n = v1.size();
    for (unsigned int i=0; i<n; i++) {
        v1[i] *= v2[i];
    }
    return v1;
}

dvector operator*(const dvector& v1, const dvector& v2)
{
    dvector v(v1);
    return (v *= v2);
}

dvector operator*(const dvector& v1, const double s)
{
    dvector v(v1);
    return (v *= s);
}

dvector operator*(const double s, const dvector& v1)
{
    dvector v(v1);
    return (v *= s);
}

dvector& operator/=(dvector& v1, const dvector& v2)
{
    if (v1.size() != v2.size()) {
        throw debugException("mathUtils::operator/=: cannot divide vectors of different sizes.");
    }
    unsigned int n = v1.size();
    for (unsigned int i=0; i<n; i++) {
        v1[i] /= v2[i];
    }
    return v1;
}

dvector operator/(const dvector& v1, const dvector& v2)
{
    dvector v(v1);
    return (v /= v2);
}

dvector& operator/=(dvector& v1, const double s)
{
    unsigned int n = v1.size();
    for (unsigned int i=0; i<n; i++) {
        v1[i] /= s;
    }
    return v1;
}

dvector operator/(const dvector& v1, const double s)
{
    dvector v(v1);
    return (v /= s);
}

vector<bool> operator>(const dvector& v1, const dvector& v2)
{
    if (v1.size() != v2.size())
        throw debugException("mathUtils::operator>: cannot compare vectors of different sizes.");

    dvector::size_type n = v1.size();
    vector<bool> v(n);
    for (dvector::size_type i=0; i<n; i++) {
        v[i] = v1[i] > v2[i];
    }
    return v;
}

vector<bool> operator<(const dvector& v1, const dvector& v2)
{
    if (v1.size() != v2.size())
        throw debugException("mathUtils::operator<: cannot compare vectors of different sizes.");

    dvector::size_type n = v1.size();
    vector<bool> v(n);
    for (dvector::size_type i=0; i<n; i++) {
        v[i] = v1[i] < v2[i];
    }
    return v;
}

vector<bool> operator>=(const dvector& v1, const dvector& v2)
{
    if (v1.size() != v2.size())
        throw debugException("mathUtils::operator>=: cannot compare vectors of different sizes.");

    dvector::size_type n = v1.size();
    vector<bool> v(n);
    for (dvector::size_type i=0; i<n; i++) {
        v[i] = v1[i] >= v2[i];
    }
    return v;
}

vector<bool> operator<=(const dvector& v1, const dvector& v2)
{
    if (v1.size() != v2.size())
        throw debugException("mathUtils::operator<=: cannot compare vectors of different sizes.");

    dvector::size_type n = v1.size();
    vector<bool> v(n);
    for (dvector::size_type i=0; i<n; i++) {
        v[i] = v1[i] <= v2[i];
    }
    return v;
}

vector<bool> operator==(const dvector& v1, const dvector& v2)
{
    if (v1.size() != v2.size())
        throw debugException("mathUtils::operator==: cannot compare vectors of different sizes.");
    dvector::size_type n = v1.size();
    vector<bool> v(n);
    for (dvector::size_type i=0; i<n; i++) {
        v[i] = v1[i] == v2[i];
    }
    return v;
}

vector<bool> operator!=(const dvector& v1, const dvector& v2)
{
    if (v1.size() != v2.size())
        throw debugException("mathUtils::operator!=: cannot compare vectors of different sizes.");

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
    if (n != v2.size())
        throw debugException("mathUtils::operator&&: cannot compare vectors of different sizes.");

    vector<bool> out(n);
    for (vector<bool>::size_type i=0; i<n; i++) {
        out[i] = v1[i] && v2[i];
    }
    return out;
}

vector<bool> operator||(const vector<bool>& v1, const vector<bool>& v2)
{
    vector<bool>::size_type n = v1.size();
    if (n != v2.size())
        throw debugException("mathUtils::operator||: cannot compare vectors of different sizes.");

    vector<bool> out(n);
    for (vector<bool>::size_type i=0; i<n; i++) {
        out[i] = v1[i] || v2[i];
    }
    return out;
}
