#pragma once
#include <iostream>
#include <chrono>
#include <thread>
using namespace std;
inline int uwait(double usec)
{
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    while (chrono::duration<double, std::micro>(end - start).count() < usec) {
        end = std::chrono::high_resolution_clock::now();
    }
    return (end - start).count();
}
class Timer {
public:
    chrono::steady_clock::time_point start;
    chrono::steady_clock::time_point loop;
    chrono::steady_clock::time_point end;
    chrono::duration<double, std::micro> interval;
    Timer() {
        start = std::chrono::high_resolution_clock::now();
        loop = std::chrono::high_resolution_clock::now();
        end = std::chrono::high_resolution_clock::now();
        interval = (loop - start);
    }
    ~Timer() {}
    void set_to_zero() {
        start = loop = end = std::chrono::high_resolution_clock::now();
        interval = (loop - start);
    }
    double get_loop_interval() {
        end = std::chrono::high_resolution_clock::now();
        interval = (end - loop);
        loop = end;
        return interval.count();
    }
    double get_full_interval() {
        end = std::chrono::high_resolution_clock::now();
        interval = (end - start);
        return interval.count();
    }
}; 
