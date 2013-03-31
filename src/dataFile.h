#pragma once

// Force use of the 1.8 HDF5 API
#define H5Dcreate_vers 2
#define H5Dopen_vers 2

#include <string>
#include <hdf5.h>

#include "mathUtils.h"

//! Represents an HDF5 data file
class DataFile
{
public:
    DataFile();

    //! Create and open an HDF5 file.
    //! @copydetails open()
    explicit DataFile(const std::string& filename,
                      std::ios_base::openmode mode=std::ios_base::app);

    ~DataFile();

    //! Add a single scalar named *name* with value *x*.
    void writeScalar(const std::string& name, const double x);

    //! Read a scalar named *name*
    double readScalar(const std::string& name);

    //! Write a 1D array (`std::vector<double>`) named *name* with values *v*.
    void writeVector(const std::string& name, const dvector& v);

    //! Read the contents of a 1D array named *name* as a `std::vector<double>`.
    dvector readVector(const std::string& name);

    //! Write a 1D array (`Eigen::ArrayXd`) named *name* with values *v*.
    void writeVec(const std::string& name, const dvec& v);

    //! Read the contents of a 1D array named *name* as an `Eigen::ArrayXd`.
    dvec readVec(const std::string& name);

    //! Write a 2D array (`Eigen::ArrayXXd`) named *name* with values *y*.
    //! Set `transpose` to `true` to change the orientation of the matrix.
    void writeArray2D(const std::string& name, const dmatrix& y, bool transpose=false);

    //! Read a 2D array (`Eigen::ArrayXXd`) named *name*.
    //! Set `transpose` to `true` to change the orientation of the matrix.
    dmatrix readArray2D(const std::string& name, bool transpose=false);

    //! Open an HDF5 data file.
    //! @param filename The name of the file. The file will be created if it
    //!     does not already exist.
    //! @param mode May be either `std::ios_base::app` to append to an
    //!     existing file or `std::ios_base::trunc` to truncate and discard
    //!     any existing contents.
    void open(const std::string& filename,
              std::ios_base::openmode mode=std::ios_base::app);

    //! Close the associated file.
    void close();

private:
    void require_file_open();
    hid_t get_dataset(const std::string& name, hid_t datatype, hid_t dataspace);

    std::string filename;
    hid_t file;
    bool fileIsOpen;
};

void dataFileTest();
