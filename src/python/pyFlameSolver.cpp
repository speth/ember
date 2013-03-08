#include <boost/python.hpp>

#include "../flameSolver.h"
#include "boost/numpy/eigen.hpp"

using namespace boost::python;

void exportFlameSolver()
{
    class_<FlameSolver, boost::noncopyable>("FlameSolver")
        .def("setOptions", &FlameSolver::setOptions)
        .def("initialize", &FlameSolver::initialize)
        .def("step", &FlameSolver::step)
        .def("finalize", &FlameSolver::finalize)
        .def("writeStateFile", &FlameSolver::writeStateFile)
        .def("writeTimeseriesFile", &FlameSolver::writeTimeseriesFile)
        .def_readonly("timeVector", &FlameSolver::timeVector)
        .def_readonly("heatReleaseRate", &FlameSolver::heatReleaseRate)
        .def_readonly("consumptionSpeed", &FlameSolver::consumptionSpeed)
        .def_readonly("flamePosition", &FlameSolver::flamePosition)
        .def_readonly("terminationCondition", &FlameSolver::terminationCondition)
        .def_readonly("grid", &FlameSolver::grid)
        .def_readonly("U", &FlameSolver::U)
        .def_readonly("T", &FlameSolver::T)
        .def_readonly("Y", &FlameSolver::Y)
        .def_readonly("dUdtDiff", &FlameSolver::dUdtDiff)
        .def_readonly("dUdtConv", &FlameSolver::dUdtConv)
        .def_readonly("dUdtProd", &FlameSolver::dUdtProd)
        .def_readonly("dTdtDiff", &FlameSolver::dTdtDiff)
        .def_readonly("dTdtConv", &FlameSolver::dTdtConv)
        .def_readonly("dTdtProd", &FlameSolver::dTdtProd)
        .def_readonly("dYdtDiff", &FlameSolver::dYdtDiff)
        .def_readonly("dYdtConv", &FlameSolver::dYdtConv)
        .def_readonly("dYdtProd", &FlameSolver::dYdtProd)
        .def_readonly("drhodt", &FlameSolver::drhodt)
        .def_readonly("rho", &FlameSolver::rho)
        .def_readonly("jCorr", &FlameSolver::jCorr)
        .def_readonly("sumcpj", &FlameSolver::sumcpj)
        .def_readonly("qDot", &FlameSolver::qDot)
        .def_readonly("wDot", &FlameSolver::wDot)
        .def_readonly("Wmx", &FlameSolver::Wmx)
        .def_readonly("W", &FlameSolver::W)
        .def_readonly("mu", &FlameSolver::mu)
        .def_readonly("k", &FlameSolver::lambda)
        .def_readonly("cp", &FlameSolver::cp)
        .def_readonly("cpSpec", &FlameSolver::cpSpec)
        .def_readonly("rhoD", &FlameSolver::rhoD)
        .def_readonly("Dkt", &FlameSolver::Dkt)
        .def_readonly("hk", &FlameSolver::hk)
        .def_readonly("jFick", &FlameSolver::jFick)
        .def_readonly("jSoret", &FlameSolver::jSoret)
        .def_readonly("dYdtCross", &FlameSolver::dTdtCross)
        ;

    class_<OneDimGrid, boost::noncopyable>("OneDimGrid")
        .def_readonly("x", &OneDimGrid::x)
        ;
}
