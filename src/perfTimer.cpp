#include "perfTimer.h"

// Platform-independent functions

unsigned long PerfTimer::getCallCount() const
{
    return callCount;
}

#ifdef _WIN32 // Windows-specific functions

PerfTimer::PerfTimer()
{
    callCount = 0;
    cummulativeTime = 0;
    LARGE_INTEGER f;
    QueryPerformanceFrequency(&f);
    clockFreq = double(f.QuadPart);
}

double PerfTimer::getTime() const
{
    return double(cummulativeTime)/clockFreq;
}

void PerfTimer::reset()
{
    callCount = 0;
    cummulativeTime = 0;
}

void PerfTimer::start()
{
    callCount++;
    QueryPerformanceCounter(&t1);
}

void PerfTimer::stop()
{
    QueryPerformanceCounter(&t2);
    cummulativeTime += t2.QuadPart - t1.QuadPart;
}

void PerfTimer::resume()
{
    QueryPerformanceCounter(&t1);
}

#else // POSIX implementation

PerfTimer::PerfTimer()
{
    callCount = 0;
    cummulativeTime = 0;
}

double PerfTimer::getTime() const
{
    return ((double) cummulativeTime)/1e6;
}

void PerfTimer::reset()
{
    callCount = 0;
    cummulativeTime = 0;
}

void PerfTimer::start()
{
    callCount++;
    gettimeofday(&t1,NULL);
}

void PerfTimer::stop()
{
    gettimeofday(&t2,NULL);
    cummulativeTime += (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
}

void PerfTimer::resume()
{
    gettimeofday(&t1,NULL);
}

#endif
