#include "perfTimer.h"

// Platform-independent functions

namespace {

unsigned long int zeroInit() {
    return 0;
}

template <class T>
T accumulate(T a, T b)
{
    return a+b;
}

} // local namespace

unsigned long PerfTimer::getCallCount() const
{
    tbb::combinable<unsigned long int> tmp(callCount);
    return tmp.combine(accumulate<unsigned long int>);
}

void PerfTimer::reset()
{
    callCount.clear();
    cumulativeTime.clear();
}

void PerfTimer::start()
{
    callCount.local()++;
    resume();
}

#ifdef _WIN32 // Windows-specific functions

LARGE_INTEGER largeIntegerZero() {
    LARGE_INTEGER x;
    x.QuadPart = 0;
    return x;
}

PerfTimer::PerfTimer()
    : callCount(&zeroInit)
    , cumulativeTime(&zeroInit)
    , t1(&largeIntegerZero)
{
    LARGE_INTEGER f;
    QueryPerformanceFrequency(&f);
    clockFreq = double(f.QuadPart);
}

double PerfTimer::getTime() const
{
    tbb::combinable<LONGLONG> tmp(cumulativeTime);
    return double(tmp.combine(accumulate<LONGLONG>))/clockFreq;
}

void PerfTimer::stop()
{
    LARGE_INTEGER t2;
    QueryPerformanceCounter(&t2);
    cumulativeTime.local() += t2.QuadPart - t1.local().QuadPart;
}

void PerfTimer::resume()
{
    QueryPerformanceCounter(&t1.local());
}

#else // POSIX implementation

PerfTimer::PerfTimer()
    : callCount(&zeroInit)
    , cumulativeTime(&zeroInit)
{
}

double PerfTimer::getTime() const
{
    tbb::combinable<unsigned long int> tmp(cumulativeTime);
    return ((double) tmp.combine(accumulate<unsigned long int>))/1e6;
}

void PerfTimer::stop()
{
    timeval t2;
    timeval& t1loc = t1.local();
    gettimeofday(&t2,NULL);
    cumulativeTime.local() += (t2.tv_sec-t1loc.tv_sec)*1000000 + (t2.tv_usec-t1loc.tv_usec);
}

void PerfTimer::resume()
{
    gettimeofday(&t1.local(), NULL);
}

#endif
