#include "mathUtils.h"
#include "sundialsUtils.h"
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
        throw DebugException("mathUtils::maxval: bad range specified");

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
        throw DebugException("mathUtils::minval: bad range specified");

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
        throw DebugException("mathUtils::sum: bad range specified");

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
        throw DebugException("mathUtils::mean: bad range specified");

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

dvector mathUtils::abs(const dvector& v)
{
    dvector a(v);
    for (dvector::size_type i=0; i<v.size(); i++)
    {
        a[i] = std::abs(a[i]);
    }
    return a;
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

dvector mathUtils::linspace(const double x1, const double x2, const int n)
{
    dvector v(n);
    for (int i=0; i<n; i++) {
        v[i] = x1 + ((x2-x1)*i)/(n-1);
    }
    return v;
}


dvec mathUtils::splines(const dvec& xIn, const dvec& yIn, const dvec& xOut)
{
    if (xIn.rows() != yIn.rows()) {
        throw DebugException((boost::format(
            "mathUtils::splines: error: xIn (%i) and yIn (%i) must be the same size.") %
            xIn.rows() % yIn.rows()).str());
    }

    int nOut = xOut.rows();
    int nIn = xIn.rows();

    vector<dvec> c = computeSplines(xIn, yIn);

    // Evaluate the spline at the selected points
    double dx;
    dvec yOut(nOut);

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

double mathUtils::splines(const dvec& xIn, const dvec& yIn, const double xOut)
{
    if (xIn.rows() != yIn.rows()) {
        throw DebugException((boost::format(
            "mathUtils::splines: error: xIn (%i) and yIn (%i) must be the same size.") %
            xIn.rows() % yIn.rows()).str());
    }

    int nIn = xIn.rows();

    vector<dvec> c = computeSplines(xIn, yIn);

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
        throw DebugException("mathUtils::operator+=: cannot add vectors of different sizes.");

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
        throw DebugException("mathUtils::operator-=: cannot subtract vectors of different sizes.");

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
        throw DebugException("mathUtils::operator*=: cannot multiply vectors of different sizes.");
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
        throw DebugException("mathUtils::operator/=: cannot divide vectors of different sizes.");
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
        throw DebugException("mathUtils::operator>: cannot compare vectors of different sizes.");

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
        throw DebugException("mathUtils::operator<: cannot compare vectors of different sizes.");

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
        throw DebugException("mathUtils::operator>=: cannot compare vectors of different sizes.");

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
        throw DebugException("mathUtils::operator<=: cannot compare vectors of different sizes.");

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
        throw DebugException("mathUtils::operator==: cannot compare vectors of different sizes.");
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
        throw DebugException("mathUtils::operator!=: cannot compare vectors of different sizes.");

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
        throw DebugException("mathUtils::operator&&: cannot compare vectors of different sizes.");

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
        throw DebugException("mathUtils::operator||: cannot compare vectors of different sizes.");

    vector<bool> out(n);
    for (vector<bool>::size_type i=0; i<n; i++) {
        out[i] = v1[i] || v2[i];
    }
    return out;
}
