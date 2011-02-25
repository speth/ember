#include "gilReleaser.h"

#include <Python.h>
#include <iostream>

__thread bool gil_is_released = false;

GILReleaser::GILReleaser()
{
    if (!gil_is_released && PyEval_ThreadsInitialized()) {
        save = PyThreadState_Swap(NULL);
        PyEval_ReleaseLock();
        reacquire = true;
        gil_is_released = true;
    } else {
        reacquire = false;
    }
}

GILReleaser::~GILReleaser()
{
    if (reacquire) {
        PyEval_AcquireLock();
        PyThreadState_Swap(reinterpret_cast<PyThreadState*>(save));
        gil_is_released = false;
    }
}
