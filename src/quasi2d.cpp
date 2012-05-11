#include "quasi2d.h"

#include "dataFile.h"

void BilinearInterpolator::open(const std::string& filename,
                                const std::string& path,
                                const std::string& xcoords,
                                const std::string& ycoords)
{
    DataFile file("filename");
    data_ = file.readArray2D(path);
    x_ = file.readVector(xcoords);
    y_ = file.readVector(ycoords);

    for (size_t i = 0; i < x_.size() - 1; i++) {
        xi_[x_[i]] = i;
    }
    xi_[-1e300] = 0;
    xi_[1e300] = x_.size() - 2;

    for (size_t i = 0; i < y_.size() - 1; i++) {
        yi_[y_[i]] = i;
    }
    yi_[-1e300] = 0;
    yi_[1e300] = y_.size() - 2;
}

double BilinearInterpolator::get(double x, double y) const
{
    typedef std::map<double, size_t>::const_iterator iter_t;
    size_t i = xi_.upper_bound(x)->second;
    size_t j = yi_.upper_bound(y)->second;
    double x1 = x_[i];
    double x2 = x_[i+1];
    double y1 = y_[j];
    double y2 = y_[j+1];
    return (data_(i,j) * (x2-x)* (y2-y) +
        data_(i+1,j) * (x-x1) * (y2-y) +
        data_(i,j+1) * (x2-x) * (y-y1) +
        data_(i+1,j+1) * (x-x1) * (y-y1)) / ((x2-x1) * (y2-y1));
}
