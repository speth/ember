#include "matlabFile.h"

matFile::matFile(void) 
	 : file(NULL)
{
}

matFile::~matFile(void)
{
	close();
}

matFile::matFile(const std::string& filename)
{
	open(filename);
}

void matFile::close(void)
{
	if (file) {
		Mat_Close(file);
		file = NULL;
	}
}

bool matFile::open(const std::string& filename)
{
	if (boost::filesystem::exists(filename)) {
		file = Mat_Open(filename.c_str(), MAT_ACC_RDWR);
	} else {
		file = Mat_Create(filename.c_str(), NULL);
	}
	
	return (file != NULL) ? true : false;
}

void matFile::writeVector(const std::string& name, const dvector& v) 
{
	if (file) {

		int dims[2] = {v.size(),1};
		matvar_t* var = Mat_VarCreate(name.c_str(),MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, (void*) &v[0], 0);
		var->compression = COMPRESSION_ZLIB;
		Mat_VarWrite(file, var, COMPRESSION_ZLIB);
		Mat_VarFree(var);
	}
}

dvector matFile::readVector(const std::string& name)
{
	dvector v;
	if (!file) {
		return v;
	}

	// The fact that Mat_VarReadInfo doesn't take a const char* is really annoying
	char* varName = new char[name.size()];
	strcpy(varName,name.c_str());
	matvar_t* var = Mat_VarReadInfo(file,varName);
	int stride[2] = {1,1};
	int start[2] = {0,0};
	int edge[2] = {var->dims[0],var->dims[1]};

	size_t classSize = Mat_SizeOfClass(var->class_type);
	if (classSize != SIZEOF_DOUBLE) {
		std::cout << "matFile::readVector: Error: Data type is not double" << std::endl;
		throw;
	}

	v.resize(edge[0]*edge[1]);
	int fail = Mat_VarReadData(file, var, &v[0], start, stride, edge);
	if (fail) {
		std::cout << "matFile::readVector: Error reading variable \"" << name << "\" from \"" << filename << "\"." << std::endl;
	}

	Mat_VarFree(var);
	delete &varName; // is this right?
	return v;
}