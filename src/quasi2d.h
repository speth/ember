#pragma once

#include "mathUtils.h"

#include <string>
#include <map>

//! Bilinear interpolation from data on a rectangular mesh.
//! The mesh spacing can be irregular in each direction.
class BilinearInterpolator
{
public:
    //! Read the data from a variable saved in an HDF5 file
    void open(const std::string& filename, const std::string& path,
              const std::string& xcoords="x",
              const std::string& ycoords="y");

    //! Interpolate at point x,y
    double get(double x, double y) const;

private:
    dmatrix data_; //!< Data values on a rectangular grid
    dvector x_, y_; //!< x and y values defining the grid

    //! Indices corresponding to x and y values to speed up interpolation
    std::map<double, size_t> xi_, yi_;
};
