#pragma once
#include <sys/time.h>
#include <time.h>

class perfTimer {
private:
    unsigned long int callCount;
    unsigned long int cummulativeTime; // stored in microseconds

    timeval t1, t2;

public:
    perfTimer(void);
    ~perfTimer(void) {}

    void start(void);
    void stop(void);
    void reset(void);
    void resume(void); // starts timer without incrementing callCount
    double getTime(void) const; // elapsed time in seconds
    unsigned long int getCallCount(void) const;
};

