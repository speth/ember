#include <boost/python.hpp>

#include "../src/flameSolver.h"

using namespace boost::python;

void exportFlameSolver()
{
    class_<FlameSolver>("FlameSolver")
        .def(init<object>())
        .def("initialize", &FlameSolver::initialize)
        .def("run", &FlameSolver::tryrun)
        .def("writeStateFile", &FlameSolver::writeStateFile)
        .def_readonly("timeVector", &FlameSolver::timeVector)
        .def_readonly("heatReleaseRate", &FlameSolver::heatReleaseRate)
        .def_readonly("consumptionSpeed", &FlameSolver::consumptionSpeed)
        .def_readonly("flamePosition", &FlameSolver::flamePosition)
        ;
}
