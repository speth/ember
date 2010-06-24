#pragma once

#include "mathUtils.h"
#include "sundialsUtils.h"
#include "chemistry0d.h"
#include "grid.h"

class DiffusionSystem : public sdODE, public GridBased
{
    // This is a system representing diffusion of a single solution component
public:
    // Evaluate an ODE of the form:
    // ydot = B * d/dx(D * dy/dx) + constantTerm
    int f(realtype t, sdVector& y, sdVector& ydot);

    int bandedJacobian(realtype t, sdVector& y, sdVector& ydot, sdBandMatrix& J);

    // The current solution vector
    dvector y;

    // the coefficients of the ODE
    dvector B;
    dvector D;

private:
    // The sum of the terms held constant for this system
    dvector constantTerm;
};
