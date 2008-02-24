#pragma once
#include <matio.h>
#include <string>
#include "boost/filesystem.hpp"
#include "mathUtils.h"

class matFile {
public:
	matFile(void);
	matFile(const std::string& filename);
	~matFile(void);

	void writeVector(const std::string& name, const dvector& v);
	dvector readVector(const std::string& name);

	bool open(const std::string& filename);
	void close(void);
private:
	std::string filename;
	mat_t* file;

};