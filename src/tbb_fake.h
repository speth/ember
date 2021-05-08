#pragma once

#include "debugUtils.h"

namespace tbb {

class global_control {
public:
    const static int max_allowed_parallelism = 0;
    global_control(int param, int nThreads) {
        if (nThreads > 1) {
            logFile.write("Warning: Running in single-threaded mode because TBB"
                " support was not enabled when compiling Ember.");
        }
    }
};

template <class T>
class combinable {
public:
    template <class U>
    combinable(U (*func)()) : val(func()) {}
    T combine(T (*func)(T, T)) { return val; }
    void clear() {}
    T& local() { return val; }
private:
    T val;
};

template <class T>
class enumerable_thread_specific {
public:
    enumerable_thread_specific() {}
    enumerable_thread_specific(T (*func)()) {}
    T& local() {
        return val;
    }
private:
    T val;
};

template <class T>
class blocked_range {
public:
    blocked_range(T start_, T stop_, T ignored) : start(start_), stop(stop_) {}
    T begin() const { return start; }
    T end() const { return stop; }
    T start;
    T stop;
};

template <class T>
void parallel_for(blocked_range<size_t> range, T func) {
    func(range);
};

}
