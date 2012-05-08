#pragma once

// Force use of the 1.8 HDF5 API
#define H5Dcreate_vers 2
#define H5Dopen_vers 2

#include <string>
#include <hdf5.h>

#include "mathUtils.h"

class DataFile
{
public:
    DataFile(void);
    DataFile(const std::string& filename);
    ~DataFile(void);

    void writeScalar(const std::string& name, const double x);
    double readScalar(const std::string& name);

    void writeVector(const std::string& name, const dvector& v);
    dvector readVector(const std::string& name);

    void writeArray2D(const std::string& name, const dmatrix& y, bool transpose=false);
    dmatrix readArray2D(const std::string& name, bool transpose=false);

    void open(const std::string& filename);
    void close();

private:
    void require_file_open();
    hid_t get_dataset(const std::string& name, hid_t datatype, hid_t dataspace);

    std::string filename;
    hid_t file;
    bool fileIsOpen;
};

void dataFileTest();
