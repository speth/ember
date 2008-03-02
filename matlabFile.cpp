#include "matlabFile.h"
#include "boost/filesystem.hpp"

matlabFile::matlabFile(void) 
	 : file(NULL)
{
}

matlabFile::~matlabFile(void)
{
	close();
}

matlabFile::matlabFile(const std::string& filename)
{
	open(filename);
}

void matlabFile::close(void)
{
	if (file) {
		matClose(file);
		file = NULL;
	}
}

bool matlabFile::open(const std::string& filename)
{
	myFilename = filename;
	if (boost::filesystem::exists(filename)) {
		file = matOpen(filename.c_str(), "u");
		accessModeIsUpdate = true;
	} else {
		file = matOpen(filename.c_str(), "wz");
		accessModeIsUpdate = false;
	}
	
	return (file != NULL) ? true : false;
}

void matlabFile::writeVector(const std::string& name, const dvector& v) 
{
	if (file) {
		mxArray* var;
		var = mxCreateDoubleMatrix(1,v.size(),mxREAL);
		if (v.size() != 0) {
			memcpy( mxGetPr(var), &v[0], v.size()*sizeof(double));
		}
		int status = matPutVariable(file, name.c_str(), var);
		if (status) {
			std::cout << "matFile::writeVector: Error writing variable \"" << 
				name << "\" to \"" << myFilename << "\"." << std::endl;
		}
		mxDestroyArray(var);
	}
}

dvector matlabFile::readVector(const std::string& name)
{
	if (!accessModeIsUpdate) {
		reopenForUpdate();
	}

	dvector v;
	if (!file) {
		return v;
	}
	mxArray* var = matGetVariable(file, name.c_str());
	if (var == NULL) {
		std::cout << "matFile::readVector: Error reading variable \"" << name << "\" from \"" << myFilename << "\"." << std::endl;
	}

	if (!mxIsDouble(var)) {
		std::cout << "matFile::readVector: Error: Data type is not double" << std::endl;
		throw;
	}

	int n = mxGetNumberOfElements(var);
	v.resize(n);
	memcpy(&v[0], mxGetPr(var), n*sizeof(double));

	mxDestroyArray(var);

	return v;
}

void matlabFile::writeArray2D(const std::string& name, const Cantera::Array2D& y) 
{

	if (file) {
		mxArray* var;
		var = mxCreateDoubleMatrix(y.nRows(), y.nColumns(), mxREAL);
		if (y.data().size() != 0) {
			memcpy( mxGetPr(var), &y.data()[0], y.data().size()*sizeof(double));
		}
		int status = matPutVariable(file, name.c_str(), var);
		if (status) {
			std::cout << "matFile::writeVector: Error writing variable \"" << 
				name << "\" to \"" << myFilename << "\"." << std::endl;
		}
		mxDestroyArray(var);
	}
}

Cantera::Array2D matlabFile::readArray2D(const std::string& name)
{
	if (!accessModeIsUpdate) {
		reopenForUpdate();
	}

	Cantera::Array2D y;
	if (!file) {
		return y;
	}
	mxArray* var = matGetVariable(file, name.c_str());
	if (var == NULL) {
		std::cout << "matlabFile::readArray2D: Error reading variable \"" << name << "\" from \"" << myFilename << "\"." << std::endl;
	}

	if (!mxIsDouble(var)) {
		std::cout << "matlabFile::readArray2D: Error: Data type is not double." << std::endl;
		throw;
	}

	int nDims = mxGetNumberOfDimensions(var);
	if (nDims != 2) {
		std::cout << "matlabFile::readArray2D: Error: Matlab array is not 2D." << std::endl;
	}

	int n = mxGetN(var);
	int m = mxGetM(var);
	y.resize(n,m);
	memcpy(&y(0,0), mxGetPr(var), n*m*sizeof(double));
	mxDestroyArray(var);

	return y;
}


void matlabFile::reopenForUpdate(void)
{
	close();
	open(myFilename);
	
}