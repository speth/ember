#pragma once

//#include <Python.h>
#ifdef _MSC_VER
#define __thread __declspec(thread)
#endif

class GILReleaser {
    // To release the Python GIL in a function, create an instance of this
    // class in that function. The constructor will release the GIL, and
    // the destructor will automatically reacquire it.
public:
    GILReleaser();
    ~GILReleaser();
private:
    void* save;
    bool reacquire;
};
