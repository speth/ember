#include "module.h"
#include "pyConverters.h"
#include "boost/numpy.hpp"

using namespace boost::python;

BOOST_PYTHON_MODULE(_pyro)
{
    boost::numpy::initialize();
    exportConverters();
    exportContainers();

    exportFlameSolver();
    exportReadConfig();
}
