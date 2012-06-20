// -*- Mode: C++ -*-
#pragma once

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#else
#include <sys/time.h>
#include <time.h>
#endif

class PerfTimer {
private:
    unsigned long int callCount;

#ifdef _WIN32
    double clockFreq;
    LARGE_INTEGER t1, t2;
    LONGLONG cummulativeTime; // stored in units of queryPerformanceCounter
#else
    unsigned long int cummulativeTime; // stored in microseconds
    timeval t1, t2;
#endif

public:
    PerfTimer();
    ~PerfTimer() {};

    void start();
    void stop();
    void reset();
    void resume(); // starts timer without incrementing callCount
    double getTime() const; // elapsed time in seconds
    unsigned long int getCallCount() const;
};
