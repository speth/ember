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
        ;
}
