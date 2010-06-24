#pragma once

#include "mathUtils.h"
#include "sundialsUtils.h"
#include "chemistry0d.h"
#include "grid.h"
#include "integrator.h"

namespace BoundaryCondition {
    enum BC { FixedValue, ZeroGradient };
}

class DiffusionSystem : public LinearODE, public GridBased
{
    // This is a system representing diffusion of a single solution component,
    // represented by an ODE in one of the following forms:
    //     ydot = B * d/dx(D * dy/dx) + C
    //     ydot = B/r * d/dr(r * D * dy/dr) + C

public:

    // The Jacobian matrix assoicated with the ODE, requested by the integrator
    void get_J(sdBandMatrix& J);

    // provides the constant term to the integrator
    void get_c(dvector& y);

    // The current solution vector
    dvector y;

    // the coefficients of the ODE
    dvector B; // pre-factor
    dvector C; // constant term
    dvector D; // "diffusion" coefficient

    BoundaryCondition::BC leftBC;
    BoundaryCondition::BC rightBC;

private:
};
