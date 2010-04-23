#pragma once
#include <string>
#include "mat.h"
#include "mathUtils.h"

class matlabFile {
public:
    matlabFile(void);
    matlabFile(const std::string& filename);
    ~matlabFile(void);

    void writeScalar(const std::string& name, const double x);
    double readScalar(const std::string& name);

    void writeVector(const std::string& name, const dvector& v);
    dvector readVector(const std::string& name);

    void writeArray2D(const std::string& name, const Cantera::Array2D& y);
    Cantera::Array2D readArray2D(const std::string& name);



    bool open(const std::string& filename);
    void close(void);
private:
    std::string myFilename;
    MATFile* file;
    bool accessModeIsUpdate;

    void reopenForUpdate(void);

};
