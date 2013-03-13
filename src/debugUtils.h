#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <exception>

#include <boost/format.hpp>
using boost::format;

class LogFile
{
public:
    LogFile();
    explicit LogFile(const std::string& filename);
    ~LogFile();

    void open(const std::string& filename);
    void close();

    template <class T>
    void write(const T& other, bool newline=true) {
        if (haveFile) {
            file << other;
            if (newline) {
                file << std::endl;
            }
        } else {
            std::cout << other;
            if (newline) {
                std::cout << std::endl;
            }
#ifndef NDEBUG
            std::cout.flush();
            std::cerr.flush();
#endif
        }
    }

private:
    bool haveFile;
    std::ofstream file;
};

class debugParameters
{
public:
    static void setParameters(bool adapt, bool regrid, bool timesteps,
                              bool radius, bool verbose);
    static bool debugAdapt;
    static bool debugRegrid;
    static bool debugTimesteps;
    static bool debugFlameRadiusControl;
    static bool veryVerbose;
};

extern LogFile logFile;

class DebugException : public std::exception
{
public:
    std::string errorString;
    DebugException(void);
    ~DebugException(void) throw() {}
    DebugException(const std::string& error);
    virtual const char* what() const throw();
};
