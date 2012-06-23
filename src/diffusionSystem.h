#pragma once

#include "mathUtils.h"
#include "chemistry0d.h"
#include "grid.h"
#include "integrator.h"

class DiffusionSystem : public TridiagonalODE, public GridBased
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
    void get_A(dvec& a, dvec& b, dvec& c);

    // Provides the constant term to the integrator
    void get_k(dvec& k);

    void resize(size_t N);
    void resetSplitConstants();

    // the coefficients of the ODE
    dvec B; // pre-factor
    dvec D; // "diffusion" coefficient

    // Balancing constant introduced by the splitting method
    dvec splitConst; // constant terms

    // Number of points in the region to be solved (= jRight - jLeft + 1)
    size_t N;

    // Parameters for wall flux boundary condition
    double yInf;
    double wallConst;

private:
    // constants used to compute the matrix coefficients
    dvec c1;
    dvec c2;
};
