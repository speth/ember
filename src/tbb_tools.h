#pragma once

#include "config.h"

#ifdef EMBER_USE_TBB
#include "tbb_real.h"
#else
#include "tbb_fake.h"
#endif

//! Wrapper class for calling member functions in a `tbb::parallel_for`.
template<class T>
class TbbWrapper {
public:
    typedef void (T::*RangedFunc)(size_t x, size_t y);

    //! Constructor
    //! @param func_ Pointer to a member function of class T that takes the
    //!     half-open interval `[x,y)` and operates on that subset of the
    //!     `parallel_for`'s range.
    //! @param parent_ Pointer to the object whose member function will be
    //!     called. Typically `this`.
    TbbWrapper(RangedFunc func_, T* parent_) :
        func(func_),
        parent(*parent_)
    { }

    //! Called internally by the TBB library
    void operator()(const tbb::blocked_range<size_t>& r) const
    {
        (parent.*func)(r.begin(), r.end());
    }
private:
    RangedFunc func;
    T& parent;
};
