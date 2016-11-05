#include "Python.h" // Include this first to avoid Cygwin compilation issues
#include "callback.h"

// A C++ exception that holds a Python exception so that it can be re-raised
// by translate_exception()
class CallbackError : public std::exception
{
public:
    CallbackError(void* type, void* value, void* traceback) :
        m_type((PyObject*) type),
        m_value((PyObject*) value),
        m_traceback((PyObject*) traceback)
    {}
    PyObject* m_type;
    PyObject* m_value;
    PyObject* m_traceback;
};

void LoggerCallback::eval(const std::string& name, int flag) const
{
    void* err[3] = {0, 0, 0};
    m_func(name, flag, m_pyobj, err);
    if (err[0]) {
        throw CallbackError(err[0], err[1], err[2]);
    }
}

double IntegratorCallback::eval(double x, double t, double U, double T, dvec& y) const
{
    void* err[3] = {0, 0, 0};
    double z = m_func(x, t, U, T, y, m_pyobj, err);
    if (err[0]) {
        throw CallbackError(err[0], err[1], err[2]);
    }
    return z;
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
        PyErr_Restore(exn.m_type, exn.m_value, exn.m_traceback);
    } catch (const std::out_of_range& exn) {
        PyErr_SetString(PyExc_IndexError, exn.what());
    } catch (const std::exception& exn) {
        PyErr_SetString(PyExc_Exception, exn.what());
    } catch (...) {
        PyErr_SetString(PyExc_Exception, "Unknown exception");
    }
    return 0;
}
