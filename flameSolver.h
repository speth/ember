#pragma once

#include "flameSys.h"

class flameSolver
{
public:
	void setOptions(const configOptions& options);
	void run(void);
private:
	configOptions options;
};
