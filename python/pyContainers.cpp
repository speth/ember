#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <vector>

using namespace boost::python;

void exportContainers()
{
    class_<std::vector<double> >("dvector")
        .def(vector_indexing_suite<std::vector<double> >())
        ;
}
