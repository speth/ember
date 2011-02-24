#pragma once

#include <boost/python.hpp>
#include <vector>

using namespace boost::python;

typedef std::vector<double> dvector;

struct Dvector_from_PyObject
{
    Dvector_from_PyObject();
    static void* convertible(PyObject* obj_ptr);
    static void construct(PyObject* obj_ptr, converter::rvalue_from_python_stage1_data* data);
};

