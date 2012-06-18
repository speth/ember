#pragma once

#include "tbb/mutex.h"
#include <vector>

template <class T>
class ObjectPool
{
public:
    void resize(size_t N) {
        data.clear();
        for (size_t i=0; i<N; i++) {
            data.push_back(new T);
        }
        mutexes.resize(N);
    }

    T** acquire() {
        size_t i = 0;
        while (true) {
            bool acquired = mutexes[i].try_lock();
            if (acquired) {
                return &data[i];
            }
            i++;
            if (i == data.size()) {
                i = 0;
            }
        }
    }

    void release(T** obj) {
        size_t i = obj - &data[0];
        mutexes[i].unlock();
    }

    T* operator[](size_t i) {
        return data[i];
    }

private:
    std::vector<T*> data;
    std::vector<tbb::mutex> mutexes;
};
