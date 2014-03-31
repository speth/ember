#pragma once

#include "mathUtils.h"

#include <map>

//! Bilinear interpolation from data on a rectangular mesh.
//! The mesh spacing can be irregular in each direction.
class BilinearInterpolator
{
public:
    //! Use the supplied data to set up the interpolation
    void setup(const dmatrix& data, const dvec& x, const dvec& y);

    //! Interpolate at point x,y
    double get(double x, double y) const;

private:
    //! initialize internally used data structures
    void initialize();

    dmatrix data_; //!< Data values on a rectangular grid
    dvec x_, y_; //!< x and y values defining the grid

    //! Indices corresponding to x and y values to speed up interpolation
    std::map<double, size_t> xi_, yi_;
};
