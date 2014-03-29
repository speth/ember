#include "callback.h"

#include "Python.h"

// A C++ exception that holds a Python exception so that it can be re-raised
// by translate_exception()
class CallbackError : public std::exception
{
public:
    CallbackError(void* type, void* value) :
        m_type((PyObject*) type),
        m_value((PyObject*) value)
    {}
    PyObject* m_type;
    PyObject* m_value;
};

void Callback::eval(const std::string& name) const
{
    void* err[2] = {0, 0};
    m_func(name, m_pyobj, err);
    if (err[0]) {
        throw CallbackError(err[0], err[1]);
    }
}

int translate_callback_exception()
{
    try {
        if (!PyErr_Occurred()) {
            // Let the latest Python exception pass through and ignore the
            // current one.
            throw;
        }
    } catch (CallbackError& exn) {
        // Re-raise a Python exception generated in a callback
        PyErr_SetObject(exn.m_type, exn.m_value);
    } catch (const std::out_of_range& exn) {
        PyErr_SetString(PyExc_IndexError, exn.what());
    } catch (const std::exception& exn) {
        PyErr_SetString(PyExc_Exception, exn.what());
    } catch (...) {
        PyErr_SetString(PyExc_Exception, "Unknown exception");
    }
    return 0;
}
