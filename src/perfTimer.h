#pragma once

#include "tbb_tools.h"

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#else
#include <sys/time.h>
#include <time.h>
#endif

//! A high-resolution, thread-aware timer and call counter class. Uses OS-
//! specific timers to get better precision (typically microsecond or better
//! resolution) than the cross- platform `clock` function which typically has
//! a resolution on the order of 10 ms. Multiple threads may independently
//! start() and stop() the same timer. The reported elapsed time will be the
//! sum of the elapsed time from each thread.
class PerfTimer {
public:
    PerfTimer();
    ~PerfTimer() {};

    //! Start accumulating time, and increment the call counter.
    void start();

    //! Stop accumulating time.
    void stop();

    //! Reset the accumulated time.
    void reset();

    //! Starts accumulating time without incrementing callCount.
    void resume();

    //! Get the elapsed time in seconds.
    double getTime() const;

    //! Get the total call count.
    unsigned long int getCallCount() const;

private:
    tbb::combinable<unsigned long int> callCount;

#ifdef _WIN32
    double clockFreq;
    tbb::combinable<LONGLONG> cumulativeTime; // stored in units of queryPerformanceCounter
    tbb::enumerable_thread_specific<LARGE_INTEGER> t1;
#else
    tbb::combinable<unsigned long int> cumulativeTime; // stored in microseconds
    tbb::enumerable_thread_specific<timeval> t1;
#endif
};
