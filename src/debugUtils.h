#pragma once
#include <string>
#include <iostream>
#include <exception>

class debugParameters
{
public:
    static bool debugAdapt;
    static bool debugRegrid;
    static bool debugSundials;
    static bool debugJacobian;
    static bool debugCalcIC;
    static bool debugTimesteps;
    static bool debugSolverStats;
    static bool debugPerformanceStats;
    static bool debugFlameRadiusControl;
    static bool veryVerbose;
};

class debugException : public std::exception
{
public:
    std::string errorString;
    debugException(void);
    ~debugException(void) throw() {}
    debugException(const std::string error);
};
