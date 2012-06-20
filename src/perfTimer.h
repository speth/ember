#pragma once

#include "tbb/combinable.h"
#include "tbb/enumerable_thread_specific.h"

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#else
#include <sys/time.h>
#include <time.h>
#endif

class PerfTimer {
public:
    PerfTimer();
    ~PerfTimer() {};

    void start();
    void stop();
    void reset();
    void resume(); // starts timer without incrementing callCount
    double getTime() const; // elapsed time in seconds
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
