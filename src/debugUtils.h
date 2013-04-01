#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <exception>

#include <boost/format.hpp>
using boost::format;

//! A container for global flags that set the verbosity of log messages.
class debugParameters
{
public:
    //! Set values for the various debug flags. This function exists because
    //! it's easier to expose in Cython than the individual static members of
    //! this class.
    static void setParameters(bool adapt, bool regrid, bool timesteps,
                              bool radius, bool verbose);

    //! Print verbose information while doing grid adaptation (adding and
    //! removing interior grid points).
    static bool debugAdapt;

    //! Print verbose information while regridding (adding and removing
    //! boundary grid points).
    static bool debugRegrid;

    //! Print a single summary line after each global timestep.
    static bool debugTimesteps;

    //! Print information related to the feedback control used to set the
    //! location of the flame, when this capability is active.
    static bool debugFlameRadiusControl;

    //! Print lots of extra information, primarily about the progress of the
    //! individual solvers for the split terms.
    static bool veryVerbose;
};

//! Write logging information either to a file or to `stdout` if no output
//! file has been specified. Usually accessed through the global instance
//! `logFile`.
class LogFile
{
public:
    //! Create a logger that writes to `stdout`.
    LogFile();
    //! Create a logger that writes to the file named *filename*. An existing
    //! file with the same name will be replaced.
    explicit LogFile(const std::string& filename);
    ~LogFile();

    //! Direct any further output to the file named *filename*.
    void open(const std::string& filename);

    //! Close the current log file, directing any further output to `stdout`.
    void close();

    //! Write the value of *other* to the current log destination. Templated
    //! so that any object that overloads `operator<<` for `std::ostream` will
    //! be appropriately formatted.
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

    //! Write the value of *message* to the current log destination only if
    //! the flag `debugParameters::veryVerbose` is set.
    template <class T>
    void verboseWrite(const T& message, bool newline=true) {
        if (debugParameters::veryVerbose) {
            write(message, newline);
        }
    }

private:
    bool haveFile;
    std::ofstream file;
};

//! Global LogFile instance
extern LogFile logFile;

//! Class for exceptions raised by Ember
class DebugException : public std::exception
{
public:
    std::string errorString;
    DebugException(void);
    ~DebugException(void) throw() {}
    DebugException(const std::string& error);
    virtual const char* what() const throw();
};
