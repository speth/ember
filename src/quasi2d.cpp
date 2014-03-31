#include "quasi2d.h"

void BilinearInterpolator::setup(const dmatrix& data,
                                 const dvec& x, const dvec& y)
{
    data_ = data;
    x_ = x;
    y_ = y;
    initialize();
}

void BilinearInterpolator::initialize()
{
    xi_.clear();
    yi_.clear();

    assert(data_.rows() == x_.size());
    assert(data_.cols() == y_.size());

    for (index_t i = 1; i < x_.rows() - 1; i++) {
        xi_[x_[i]] = i;
    }
    xi_[-1e300] = 0;
    xi_[1e300] = x_.rows() - 1;

    for (index_t i = 1; i < y_.rows() - 1; i++) {
        yi_[y_[i]] = i;
    }
    yi_[-1e300] = 0;
    yi_[1e300] = y_.rows() - 1;
}

double BilinearInterpolator::get(double x, double y) const
{
    size_t i = xi_.upper_bound(x)->second;
    size_t j = yi_.upper_bound(y)->second;
    double x1 = x_[i-1];
    double x2 = x_[i];
    double y1 = y_[j-1];
    double y2 = y_[j];

    return (data_(i-1,j-1) * (x2-x)* (y2-y) +
        data_(i,j-1) * (x-x1) * (y2-y) +
        data_(i-1,j) * (x2-x) * (y-y1) +
        data_(i,j) * (x-x1) * (y-y1)) / ((x2-x1) * (y2-y1));
}
