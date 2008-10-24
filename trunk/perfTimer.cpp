#include "perfTimer.h"

perfTimer::perfTimer(void)
{
	callCount = 0;
	cummulativeTime = 0;
}

perfTimer::~perfTimer(void)
{
}

void perfTimer::start(void)
{
	callCount++;
	gettimeofday(&t1,NULL);
}

void perfTimer::stop(void)
{
	gettimeofday(&t2,NULL);

	cummulativeTime += (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
}

void perfTimer::resume(void)
{
	gettimeofday(&t1,NULL);
}

void perfTimer::reset(void)
{
	callCount = 0;
	cummulativeTime = 0;
}

double perfTimer::getTime(void) const
{
	return ((double) cummulativeTime)/1e6;
}

unsigned long perfTimer::getCallCount(void) const
{
	return callCount;
}
