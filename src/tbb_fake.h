#pragma once

#include "debugUtils.h"

namespace tbb {

class task_scheduler_init {
public:
    const static int deferred = 0;
    task_scheduler_init(int ignored) {}
    void initialize(int nThreads) {
        if (nThreads > 1) {
            logFile.write("Warning: Running in single-threaded mode because TBB"
                " support was not enabled when compiling Ember.");
        }
    }
};

class mutex {
public:
    class scoped_lock {
    public:
        scoped_lock(mutex m) {};
    };
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
