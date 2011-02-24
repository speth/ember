#include <boost/python.hpp>

#include "../src/readConfig.h"

using namespace boost::python;

void exportReadConfig()
{
    class_<configOptions>("ConfigOptions")
        .def_readwrite("inputDir", &configOptions::inputDir)
        .def_readwrite("outputDir", &configOptions::outputDir)
        .def_readwrite("restartFile", &configOptions::restartFile)
        .def_readwrite("useRelativeRestartPath", &configOptions::useRelativeRestartPath)
        .def_readwrite("overrideTu", &configOptions::overrideTu)
        .def_readwrite("strainRateList", &configOptions::strainRateList)
        ;
}
