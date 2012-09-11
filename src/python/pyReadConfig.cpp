#include <boost/python.hpp>

#include "../readConfig.h"

using namespace boost::python;

void exportReadConfig()
{
    class_<ConfigOptions>("_ConfigOptions")
        .def(init<object>())
        .def_readwrite("inputDir", &ConfigOptions::inputDir)
        .def_readwrite("outputDir", &ConfigOptions::outputDir)
        .def_readwrite("restartFile", &ConfigOptions::restartFile)
        .def_readwrite("useRelativeRestartPath", &ConfigOptions::useRelativeRestartPath)
        .def_readwrite("overrideTu", &ConfigOptions::overrideTu)
        ;
}
