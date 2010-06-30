#pragma once

#include "readConfig.h"
#include "grid.h"
#include <string>
#include <cantera/kernel/Array.h>

class ProfileGenerator : public GridBased
{
    // This class creates initial flame profiles based on either parameters
    // specified in an input file, loads an existing flame profile from a
    // file, or uses combination of the two to modify the boundary conditions
    // on a profile loaded from a file.
public:
    void setOptions(configOptions& options);
    void setSourceFile(std::string filename);

    void generateProfile();
    void loadProfile();

    // Initial Profiles
    dvector U; // velocity
    dvector T; // temperature
    Cantera::Array2D Y; // mass fractions

    double t; // initial time

    // Boundary values
    double rhou, rhob, rhoLeft, rhoRight;
    double Tu, Tb, Tleft, Tright;
    dvector Yu, Yb, Yleft, Yright;

private:
    configOptions options;
};
