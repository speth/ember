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
    ProfileGenerator(configOptions& _options);
    void setSourceFile(std::string filename);

    void generateProfile();
    void loadProfile();

    // Calculate the mole fraction vector of the reactants based on the
    // equivalence ratio and the fuel and oxidizer compositions
    dvector calculateReactantMixture(void);

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
    configOptions& options;
};
