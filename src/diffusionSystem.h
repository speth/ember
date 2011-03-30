#pragma once

#include "mathUtils.h"
#include "sundialsUtils.h"
#include "chemistry0d.h"
#include "grid.h"
#include "integrator.h"

class DiffusionSystem : public LinearODE, public GridBased
{
    // This is a system representing diffusion of a single solution component,
    // represented by an ODE in one of the following forms:
    //     ydot = B * d/dx(D * dy/dx) + C
    //     ydot = B/r * d/dr(r * D * dy/dr) + C
    // The ODE in this form may be written as a linear system:
    //     ydot = Ay + C
    // where the entries of the matrix A are determined by the prefactor B,
    // the diffusion coefficients D, and the finite difference formula used.

public:
    DiffusionSystem();

    // Provide the matrix associated with the ODE to the integrator
    void get_A(sdBandMatrix& J);

    // Provides the constant term to the integrator
    void get_C(dvector& y);

    void setGrid(const oneDimGrid& grid);
    void resize(int N);
    void resetSplitConstants();

    // the coefficients of the ODE
    dvector B; // pre-factor
    dvector D; // "diffusion" coefficient

    // Diagonalized, linear approximations for terms neglected by splitting
    dvector splitConst; // constant terms
    dvector splitLinear; // diagonal jacobian components

    // Integrate only over the range [jLeft,jRight]
    // All inputs and outputs, except for the grid, cover only this range
    size_t jLeft;
    size_t jRight;

    // Number of points in the region to be solved (= jRight - jLeft + 1)
    int N;

    // Parameters for wall flux boundary condition
    double yInf;
    double wallConst;

};
