#include "debugUtils.h"

bool debugParameters::debugAdapt;
bool debugParameters::debugRegrid;
bool debugParameters::debugTimesteps;
bool debugParameters::debugFlameRadiusControl;
bool debugParameters::veryVerbose;
LogFile logFile;

LogFile::LogFile()
{
    std::cout << std::boolalpha; // prints "true" and "false" rather than 1 and 0
    haveFile = false;
}

LogFile::LogFile(const std::string& filename)
{
    open(filename);
}

LogFile::~LogFile()
{
    close();
}

void LogFile::open(const std::string& filename)
{
    if (haveFile) {
        close();
    }
    haveFile = true;
    file.open(filename.c_str());
    file << std::boolalpha; // prints "true" and "false" rather than 1 and 0
}

void LogFile::close()
{
    if (haveFile) {
        haveFile = false;
        file.close();
    }
}

debugException::debugException(void)
{
    errorString = "debugException: unspecified error.";
}

debugException::debugException(const std::string& error)
{
    errorString = error;
}

const char* debugException::what() const throw() {
    return errorString.c_str();
}
